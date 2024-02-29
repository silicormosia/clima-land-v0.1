#
#
# This file is meant for CliMA Land v0.1 example
#
#
using DataFrames: DataFrame, DataFrameRow
using Dates: isleapyear
using JLD2: load

using NetcdfIO: read_nc, save_nc!
using PkgUtility: month_days, nanmean

using Land.CanopyLayers: EVI, FourBandsFittingHybrid, NDVI, NIRv, SIF_WL, SIF_740, fit_soil_mat!
using Land.Photosynthesis: C3CLM, use_clm_td!
using Land.PlantHydraulics: VanGenuchten, create_tree
using Land.SoilPlantAirContinuum: CNPP, GPP, PPAR, SPACMono, T_VEG, initialize_spac_canopy!, prescribe_air!, prescribe_swc!, prescribe_t_leaf!, spac_beta_max, update_Cab!, update_LAI!, update_VJRWW!,
      update_par!, update_sif!, zenith_angle
using Land.StomataModels: BetaGLinearPsoil, ESMMedlyn, GswDrive, gas_exchange!, gsw_control!, prognostic_gsw!


DF_VARIABLES  = ["F_H2O", "F_CO2", "F_GPP", "SIF683", "SIF740", "SIF757", "SIF771", "NDVI", "EVI", "NIRv"];


"""

    prepare_wd(dict::Dict, wd_file::String)

Prepare weather driver dataframe to feed CliMA Land, given
- `dict` Dictionary that store grid information
- `wd_file` Weather driver file

"""
function prepare_wd(dict::Dict, wd_file::String)
    df_in = read_nc(wd_file);

    # compute T_MEAN based on the weather driver
    df_in[!,"CO2"]         .= 0.0;
    df_in[!,"Chlorophyll"] .= 0.0;
    df_in[!,"LAI"]         .= 0.0;
    df_in[!,"Vcmax"]       .= 0.0;
    df_in[!,"T_MEAN"]      .= 0.0;
    for i in eachindex(df_in.T_MEAN)
        if i < 240
            df_in[i,"T_MEAN"] = nanmean( max.(df_in.T_AIR[1:i], df_in.T_LEAF[1:i]) );
        else
            df_in[i,"T_MEAN"] = nanmean( max.(df_in.T_AIR[i-239:i], df_in.T_LEAF[i-239:i]) );
        end;
    end;

    #
    # extropolate the time series based on input variable dimensions, dimensions must be within supported settings
    #     1. extropolate the data to 1D resolution
    #     2. extropolate the data to 1H resolution
    #
    year = dict["year"];
    days = isleapyear(year) ? 366 : 365;
    @inline nt_to_1h(label::String) = (
        dat_in = dict[label];
        @assert length(dat_in) in [366, 365, 53, 52, 46, 12, 1] "Dataset length not supported";

        if length(dat_in) == 1
            dat_1d = repeat([dat_in;]; inner = days);
        elseif length(dat_in) == 12
            dat_1d = [([repeat(dat_in[_m:_m], month_days(year, _m)) for _m in 1:12]...)...]
        elseif length(dat_in) == 46
            dat_1d = repeat(dat_in; inner = 8)[1:days]
        elseif length(dat_in) in [52,53]
            dat_1d = repeat([dat_in;dat_in[end]]; inner = 7)[1:days]
        elseif length(dat_in) in [365,366]
            dat_1d = [dat_in;dat_in[end]][1:days]
        end;

        return repeat(dat_1d; inner = 24)
    );
    df_in[!,"CO2"]         .= nt_to_1h("co2_concentration");
    df_in[!,"Chlorophyll"] .= nt_to_1h("chlorophyll");
    df_in[!,"CI"]          .= nt_to_1h("clumping_index");
    df_in[!,"LAI"]         .= nt_to_1h("leaf_area_index");
    df_in[!,"Vcmax"]       .= nt_to_1h("vcmax");

    # add the fields to store outputs
    for label in DF_VARIABLES
        df_in[!,label] .= 0.0;
    end;

    return df_in
end;



"""

    prepare_spac(dict::Dict; FT = Float64)

Create a SPAC, given
- `dict` Dictionary of GriddingMachine data in a grid
- `FT` Floating number type

"""
function prepare_spac(dict::Dict; FT = Float64)
    # read general information from dict
    sm = ESMMedlyn{FT}();

    # use JULES soil depth 0.00 -- 0.10 -- 0.35 -- 1.00 -- 3.00 m, and assume 2 m deep root (z_root = -2) for all the sites
    soil_bounds = FT[0, -0.1, -0.35, -1, -3];
    z_canopy    = max(FT(0.1), dict["canopy_height"]);
    Δz          = z_canopy / 20;
    air_bounds  = collect(0:Δz:z_canopy+2*Δz);
    plant_hs    = create_tree(FT(-2), z_canopy/2, z_canopy, soil_bounds, air_bounds);

    # create a SPACMono struct, redefine the wavelength limits for PAR if ePAR is true
    node = SPACMono{FT}(soil_bounds=soil_bounds, air_bounds=air_bounds, z_canopy=z_canopy, z_root=-2, plant_hs=plant_hs, latitude=dict["latitude"], longitude=dict["longitude"], stomata_model=sm);

    for iPS in node.plant_ps
        iPS.g_min   = eps(FT);
        iPS.g_min25 = eps(FT);
        iPS.g_max   = 0.8;
        iPS.g_max25 = 0.8;
    end;

    # update soil type information per layer
    for i in eachindex(node.plant_hs.roots)
        α  = dict["soil_vg_α"][i];
        n  = dict["soil_vg_n"][i];
        Θr = dict["soil_vg_Θr"][i];
        Θs = dict["soil_vg_Θs"][i];
        node.plant_hs.roots[i].sh = VanGenuchten{FT}(stype = "JULES", α = α, n = n, Θs = Θs, Θr = Θr);
    end;

    # update leaf mass per area (from m² kg⁻¹ to g cm⁻²)
    for leaf in node.leaves_rt
        leaf.Cm = dict["leaf_mass_per_area"];
    end;

    # set up empirical model
    if typeof(sm) <: ESMMedlyn
        node.photo_set = C3CLM(FT);
        node.stomata_model.g1 = dict["g1_medlyn_c3"];
        node.stomata_model.g0 = 1e-3;
    else
        @warn "Stomatal model parameters are not initialized for $(typeof(sm))";
    end;

    # update soil color class from CLM dataset
    node.soil_opt.color = dict["soil_color"];

    # update the Vcmax, Jmax, and Vpmax
    update_VJRWW!(node, nanmean(dict["vcmax"]));

    # initialize the canopy RT model
    initialize_spac_canopy!(node);

    return node
end;


"""

Structure that store memory information

"""
Base.@kwdef mutable struct SPACMemory{FT<:AbstractFloat}
    chl::FT = -9999
    lai::FT = -9999
    vcm::FT = -9999
end;


"""

    prescribe_parameters!(spac::SPACMono{FT}, dfr::DataFrame, mem::SPACMemory{FT}, deepcopies::Vector) where {FT<:AbstractFloat}

Prescibe parameters for the SPAC, given
- `spac` Soil plant air continuum struct
- `dfr` Weather driver dataframe row
- `mem` Memory cache struct
- `deepcopies` Deepcopies of radiation used to scale direct and diffuse radiation

"""
function prescribe_parameters!(spac::SPACMono{FT}, dfr::DataFrameRow, mem::SPACMemory{FT}, deepcopies::Vector) where {FT<:AbstractFloat}
    # read the data out of dataframe row to reduce memory allocation
    df_atm::FT = dfr.P_ATM;
    df_chl::FT = dfr.Chlorophyll;
    df_cli::FT = dfr.CI;
    df_co2::FT = dfr.CO2;
    df_dif::FT = dfr.RAD_DIF;
    df_dir::FT = dfr.RAD_DIR;
    df_doy::FT = dfr.FDOY;
    df_lai::FT = dfr.LAI;
    df_sw1::FT = dfr.SWC_1;
    df_sw2::FT = dfr.SWC_2;
    df_sw3::FT = dfr.SWC_3;
    df_sw4::FT = dfr.SWC_4;
    df_tar::FT = dfr.T_AIR;
    df_tlf::FT = dfr.T_LEAF;
    df_tmn::FT = dfr.T_MEAN;
    df_vcm::FT = dfr.Vcmax;
    df_vpd::FT = dfr.VPD;
    df_wnd::FT = dfr.WIND;

    # adjust optimum t based on 10 day moving average skin temperature
    use_clm_td!(spac.photo_set, df_tmn);

    # if total LAI, Vcmax, or Chl changes, update them (add vertical Vcmax profile as well)
    trigger_lai::Bool = !isnan(df_lai) && (df_lai != mem.lai);
    trigger_vcm::Bool = !isnan(df_vcm) && (df_vcm != mem.vcm);
    trigger_chl::Bool = !isnan(df_chl) && (df_chl != mem.chl);
    if trigger_lai
        update_LAI!(spac, df_lai);
        mem.lai = df_lai;
    end;

    if trigger_lai || trigger_vcm
        update_VJRWW!(spac, df_vcm; expo = FT(0.3));
        mem.vcm = df_vcm;
    end;

    if trigger_chl
        update_Cab!(spac, df_chl; cab_2_car = FT(1/7));
        mem.chl = df_chl;
    end;

    # update clumping index
    spac.canopy_rt.Ω = df_cli;
    spac.canopy_rt.clump_a = df_cli;

    # sync the environmental conditions per layer
    prescribe_air!(spac, df_co2, df_atm, df_tar, df_vpd, df_wnd);
    prescribe_t_leaf!(spac, max(df_tar, df_tlf));

    # run the chunks below only when total radiation is higher than 10
    if df_dir + df_dif < 10
        return nothing
    end;

    # update soil water matrices per layer
    prescribe_swc!(spac, df_sw1, df_sw2, df_sw3, df_sw4);

    # update soil albedo using FourBandsFittingHybrid
    fit_soil_mat!(spac.soil_opt, spac.wl_set, spac.swc[1], FourBandsFittingHybrid());

    # update PAR related information
    spac.in_rad.E_direct  .= deepcopies[1].E_direct  .* df_dir ./ deepcopies[2];
    spac.in_rad.E_diffuse .= deepcopies[1].E_diffuse .* df_dif ./ deepcopies[3];
    spac.angles.sza = min(88, zenith_angle(spac.latitude, df_doy));
    update_par!(spac);

    return nothing
end;


"""

    run_time_step!(spac::SPACMono{FT}, dfr::DataFrame) where {FT<:AbstractFloat}

Run CliMA Land in a time step, given
- `spac` Soil plant air continuum struct
- `dfr` Weather driver dataframe row
- `ind` Time index

"""
function run_time_step!(spac::SPACMono{FT}, dfr::DataFrameRow, beta::BetaGLinearPsoil{FT}) where {FT<:AbstractFloat}
    # read the data out of dataframe row to reduce memory allocation
    df_dif::FT = dfr.RAD_DIF;
    df_dir::FT = dfr.RAD_DIR;

    # compute beta factor (based on Psoil, so canopy does not matter)
    # note here that
    #     - if the beta is applied to g1 use it in the prognostic_gsw! function
    #     - if the beta is applied to Vcmax, rescale Vcmax and use a beta = 1 in the stomatal_conductance function
    # in this example, we assumed it is applied to g1
    βm = spac_beta_max(spac, beta);

    # calculate leaf level flux per canopy layer
    for i in 1:spac.n_canopy
        iEN = spac.envirs[i];
        iPS = spac.plant_ps[i];

        # set gsw to 0 or iterate for 30 times to find steady state solution
        if df_dir + df_dif < 10
            iPS.APAR .= 0;
            iPS.g_sw .= 0;
            gsw_control!(spac.photo_set, iPS, iEN);
        else
            for _ in 1:30
                gas_exchange!(spac.photo_set, iPS, iEN, GswDrive());
                prognostic_gsw!(iPS, iEN, spac.stomata_model, βm, FT(120));
                gsw_control!(spac.photo_set, iPS, iEN);
            end;
        end;
    end;

    # calculate the SIF if there is sunlight
    if df_dir + df_dif >= 10
        update_sif!(spac);
        dfr.SIF683 = SIF_WL(spac.can_rad, spac.wl_set, FT(682.5));
        dfr.SIF740 = SIF_740(spac.can_rad, spac.wl_set);
        dfr.SIF757 = SIF_WL(spac.can_rad, spac.wl_set, FT(758.7));
        dfr.SIF771 = SIF_WL(spac.can_rad, spac.wl_set, FT(770.0));
        dfr.NDVI   = NDVI(spac.can_rad, spac.wl_set);
        dfr.EVI    = EVI(spac.can_rad, spac.wl_set);
        dfr.NIRv   = NIRv(spac.can_rad, spac.wl_set);
    end;

    # save the total flux into the DataFrame
    dfr.F_H2O = T_VEG(spac);
    dfr.F_CO2 = CNPP(spac);
    dfr.F_GPP = GPP(spac);

    return nothing
end;


"""

    run_model!(spac::SPACMono{FT}, df::DataFrame, nc_out::String) where {FT<:AbstractFloat}

Run CliMA Land at a site for the enture year, given
- `spac` Soil plant air continuum struct
- `df` Weather driver dataframe
- `nc_out` File path to save the model output

"""
function run_model!(spac::SPACMono{FT}, df::DataFrame, nc_out::String) where {FT<:AbstractFloat}
    in_rad_bak = deepcopy(spac.in_rad);
    in_dir     = in_rad_bak.E_direct' * spac.wl_set.dWL / 1000;
    in_dif     = in_rad_bak.E_diffuse' * spac.wl_set.dWL / 1000;
    deepcopies = [in_rad_bak, in_dir, in_dif];
    beta_g     = BetaGLinearPsoil{FT}();

    # set up memory
    spac_mem = SPACMemory{FT}();

    # iterate through the time steps
    for dfr in eachrow(df)
        prescribe_parameters!(spac, dfr, spac_mem, deepcopies);
        run_time_step!(spac, dfr, beta_g);
    end;

    # save simulation results to hard drive
    save_nc!(nc_out, df[:, DF_VARIABLES]);

    return nothing
end;


@time dict = load("$(@__DIR__)/debug.jld2");
@time wddf = prepare_wd(dict, "$(@__DIR__)/debug.nc");
@time spac = prepare_spac(dict);
@time run_model!(spac, wddf, "$(@__DIR__)/debug.output.nc");
