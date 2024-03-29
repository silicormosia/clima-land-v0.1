# API
```@meta
CurrentModule = Land.PlantHydraulics
```




## Plant Hydraulic System
The PlantHydraulics module provides two levels of hydraulics system:
    organ-level and plant-level. The organ-level hydraulic systems include
    Leaf, Root, and Stem (trunk and branch). The plant-level hydraulic system
    is can be any combination of the three organs (custimized definition may
    apply).




## Leaf, Root, and Stem organs
Plant hydraulics is segmented to three organ-level systems/structs
    ([`LeafHydraulics`](@ref), [`RootHydraulics`](@ref), and
    [`StemHydraulics`](@ref)) subject to an Abstract type
    ([`AbstractHydraulicOrgan`](@ref)). The major differences among the three
    structs are

- [`LeafHydraulics`](@ref) has an extra-xylary component
- [`RootHydraulics`](@ref) has a rhizosphere component
- [`RootHydraulics`](@ref) and [`StemHydraulics`](@ref) have a gravity component

See the documentation for each struct for more details:

```@docs
AbstractHydraulicOrgan
LeafHydraulics
RootHydraulics
StemHydraulics
```

To initialize a hydraulics system, one needs to provide the floating type, for
    example:

```julia
FT = Float32;
hs_leaf = LeafHydraulics{FT}();
hs_root = RootHydraulics{FT}();
hs_stem = StemHydraulics{FT}();
```




## Whole-plant organism
Plants differ in their structures, for example, some plants have a canopy far
    above the ground elevated by a trunk, some plants have a structured canopy
    supported by branch systems, and some plant has no trunk at all. To
    represent the structural differences, several types of plant hydraulics
    systems are pre-defined, and they are [`GrassLikeOrganism`](@ref),
    [`PalmLikeOrganism`](@ref), [`TreeLikeOrganism`](@ref), and
    [`TreeSimple`](@ref) structs subject to a [`AbstractPlantOrganism`](@ref)
    type. The major difference between the `HS`s are

- [`GrassLikeOrganism`](@ref) has only mutiple root and canopy layers, no trunk
    or branch
- [`PalmLikeOrganism`](@ref) has multiple root layers, a trunk, and multiple
    canopy layers, no branch system
- [`TreeLikeOrganism`](@ref) has multiple root layers, a trunk, and multiple
    branch + canopy layers, and each branch corresponds to a canopy layer
- [`TreeSimple`](@ref) has one root, one stem, and one leaf for testing purpose

See the documentation for each struct for more details:

```@docs
AbstractPlantOrganism
GrassLikeOrganism
PalmLikeOrganism
TreeLikeOrganism
TreeSimple
```

To ease the initialization of a plant hydraulics system, a few customized
    functions are provided for quick initialization. More importantly,
    modifications to each field in the struct are always allowed. The quick
    functions are [`create_grass`](@ref), [`create_palm`](@ref), and
    [`create_tree`](@ref):

```@docs
create_grass
create_palm
create_tree
```

What these functions do are to determine how many root layers and branch/canopy
    layers to add based on the tree information and environmental settings. To
    determine number of root layers, rooting depth and the soil layer
    information are required. The `z_root` is the maximal root depth in
    negative number, and `soil_bounds` is the boundaries of soil layers staring
    from 0. For example, for a `soil_bounds` of [0.0, -1.0, -2.0, -3.0, -4.0],
    a `z_root` of -1 gives 1 root layer, and a `z_root` of -1.5 or -2.0 gives 2
    root layers. The `z_trunk`, `z_canopy`, and `air_bounds` determine how many
    canopy layers to add. For example, for a `air_bounds` of [0.0, 1.0, 2.0,
    3.0, 4.0, 5.0 ... 20.0, 21.0, 22.0], a `z_trunk` of 5.0 `z_canopy` of 7.0
    give 2 canopy layers, and a `z_trunk` of 5.5 `z_canopy` of 7.0 give 2
    canopy layers. Also, the `root_index_in_soil` and `canopy_index_in_air`
    indicate which soil or air layer the root or canopy layer corresponds with,
    respectively. For instance, a index of 7 means that the canopy layer should
    use the 7th air layer.

To initialize a whole-plant hydraulic system, checkout the example below:

```julia
FT = Float32;
grass = create_grass(FT(-2.1), FT(0.5), FT(8), FT[0,-1,-2,-3], collect(FT,0:1:20));
palm  =  create_palm(FT(-2.1), FT(0.5), FT(8), FT[0,-1,-2,-3], collect(FT,0:1:20));
tree  =  create_tree(FT(-2.1), FT(0.5), FT(8), FT[0,-1,-2,-3], collect(FT,0:1:20));
treet = TreeSimple{FT}();
```




## Xylem hydraulic conductance
Plants transport water through xylem conduits (vessels in most angiosperms,
    trachieds in most gymnosperms). With the ascent of sap along the hydraulic
    system, water pressure in the conduits is typically negative. The negative
    xylem water pressure tends to pull air from surrounding tisses or the
    atmosphere into the xylem conduits, resulting in xylem cavitation. The air
    bubbles in cavitated conduits block water flow, and thus results in decline
    of water transport capability (measured by xylem hydraulic conductance).

Typically, the correlation between xylem water pressure ($P \leq 0$) and
    hydraulic conductance ($k$) is expressed by a Weibull function for
    [`WeibullSingle`](@ref) type correlation:

```math
k = k_\text{max} \cdot \exp \left( -\left( \dfrac{-P}{B} \right)^C \right)
```

where $k_\text{max}$ is the maximal hydraulic conductance, and $B$ and $C$ are
    the Weibull parameters. This correlation is also known as vulnerability
    curve (VC) to drought stress. Sometimes, plants exhibit a segmented VC, for
    example, the fibers may transport water as well and are much more resistant
    to drought than vessels. Thus, a dual Weibull function is presented for
    [`WeibullDual`](@ref) type correlation ($P \leq 0$):

```math
k = k_\text{max} \cdot \left\{ f_1 \cdot \exp \left[ -\left( \dfrac{-P}{B_1} \right)^{C_1} \right] +
                         (1 - f_1) \cdot \exp \left[ -\left( \dfrac{-P}{B_2} \right)^{C_2} \right] \right\}
```

The VC formulations are abstractized as

```@docs
AbstractXylemVC
WeibullDual
WeibullSingle
LogisticSingle
PowerSingle
```

The function to call is

```@docs
xylem_k_ratio
```

Note it here that `xylem_k_ratio(vc, p)` calculate the k without making
    temperature corrections, but the `xylem_k_ratio(vc, p_25, vis)` makes
    correction over the viscosity (the higher the viscosity, the lower the k).
    Also, `p_25` means that the pressure has been corrected to 298.15 K for
    surface tension (the higher the surface tension, the more resistant the
    xylem).

Meanwhile, there is a function to call to calculate the critical pressure,
    beyond which leaf will decicate. The critical pressure is calculated as the
    pressure at which $k$ is 0.001 of $k_\text{max}$ for
    [`WeibullSingle`](@ref) (for [`WeibullDual`](@ref), each segment need to
    reach 0.001). The functions is

```@docs
xylem_p_crit
```

Examples:

```julia
FT = Float32;
vc_1 = WeibullSingle{FT}();
vc_2 = WeibullDual{FT}();

k_1 = xylem_k_ratio(vc_1, -1.0);
k_2 = xylem_k_ratio(vc_2, -1.0);
k_3 = xylem_k_ratio(vc_1, -1.0, 1.2);
k_4 = xylem_k_ratio(vc_2, -1.0, 1.2);
```




## Rhizosphere hydraulic conductance
As mentioned above, there is a rhizosphere component in the root hydraulic
    system, and thus one needs to compute the pressure dtop along the
    rhizosphere. The soil properties are classified to [`BrooksCorey`](@ref)
    and [`VanGenuchten`](@ref) types subjected to [`AbstractSoilVC`](@ref):

```@docs
AbstractSoilVC
BrooksCorey
VanGenuchten
```

Pre-defined parameter sets are avialble, and you may quick create a soil type
    struct using [`create_soil_VC`](@ref). Note that soil type parameters are
    van Genuchten type VC, and we curve fitted the curve to provide the
    parameters for Brooks and Corey type VC.

```@docs
create_soil_VC
fit_soil_VC!
```

Correlations among soil relative water content, relative hydraulic conductance,
    and soil matrix potential are

```@docs
soil_erwc
soil_rwc
soil_swc
soil_k_ratio_erwc
soil_k_ratio_rwc
soil_k_ratio_swc
soil_k_ratio_p25
soil_p_25_erwc
soil_p_25_rwc
soil_p_25_swc
```




## Pressure and Flow
The PlantHydraulics module is designed to run numerically for the following
    reasons:

- Weibull function is cannot be integrated
- The VC is segmented, i.e., if $P > 0$, $k = k_\text{max}$ (implemented in
    [`xylem_k_ratio`](@ref))
- Once xylem cavitation occurs, it cannot be easily recovered unless $P > 0$,
    and thus there is a drought legacy effect. This is why there are a few
    fields in the [`LeafHydraulics`](@ref), [`RootHydraulics`](@ref), and
    [`StemHydraulics`](@ref) structs to store the drought history information.
- Temperature may change along the flow path. The `f_st` and `f_vis` in the
    structs help deal with these effects.

Function [`end_pressure`](@ref) calculates the xylem end pressure for an
    organ-level hysraulic system. As mentioned above, the
    [`RootHydraulics`](@ref) and [`StemHydraulics`](@ref) has a gravity
    component, and the [`RootHydraulics`](@ref) has a rhizosphere component.
    Also be aware that [`end_pressure`](@ref) accounts for temperature
    effects on surface tension and viscosity.

```@docs
end_pressure
```

Noe that function [`end_pressure`](@ref) does not update the pressure
    profiles or history in the xylem. To update these profiles, use
    [`pressure_profile!`](@ref), and to remove these legacy profiles, use
    [`inititialize_legacy!`](@ref):

```@docs
pressure_profile!
inititialize_legacy!
```

Examples:
```julia
FT = Float32;
leaf = LeafHydraulics{FT}();
p = end_pressure(leaf, FT(0.01));
@show leaf.p_element;
pressure_profile!(leaf, FT(0.01));
@show leaf.p_element;
```

## Steady state and non-steady state mode

```@docs
AbstractFlowMode
NonSteadyStateMode
SteadyStateMode
buffer_rate
```




## Root Hydraulics
Function [`end_pressure`](@ref) works for the case of only 1 root layer if
    one needs the plant base xylem water pressure. However, when there are
    multiple root layers, [`end_pressure`](@ref) does not apply. In this
    case, iterations are required to calculate the xylem end pressure for each
    root layers, and then make sure all root layers have the same xylem end
    pressure. A few functions are provided to realize this.

Function [`xylem_flow`](@ref) uses Root Solving method to calculate
    the flow rate through the [`RootHydraulics`](@ref) struct that yields the
    given xylem end pressure. The `ini` in the function is optional. However,
    using the flow rate from last instant when pressure does not differ much
    will speed up the calculation.

```@docs
xylem_flow
```

In the plant hydraulic module design, flow rate is computed for each canopy
    layer, and thus computing flow rate for each root layer is required for a
    multiple layered root system. One feasible way is to do iterations using
    [`xylem_flow`](@ref) function, i.e., iterate the xylem end
    pressure til the total flow rate equals the given value. However, this
    method is too inefficient. What I did is

- Calculate the xylem end pressure and whole root layer conductance from the
    initial flow rate;
- Calculate the mean xylem end pressure, and tune the flow rates in all root
    layers using the difference from mean pressure and root conductance;
- Sum up the new flow rates, and calculate the difference with given total flow
    rate;
- Use the calculated condutcance to weight out the differences.

The functions provided by PlantHydraulics module are

```@docs
root_pk
roots_flow!
```

However, the steps above are only 1 iteration, and can only be used for the
    non-steady state version of model. For the steady-state flow rates,
    function [`roots_flow!`](@ref) does thw work. What the function does is to
    iterate [`roots_flow!`](@ref) till the difference among the
    calculated end pressures is small enough. I also emphasize that to speed up
    the code, 3 containers are added to the [`AbstractPlantOrganism`](@ref)
    structs, and they are `cache_k`, `cache_p`, and `cache_q`.

Example:

```julia
FT = Float32;
palm = create_palm(FT(-2.1), FT(0.5), FT(8), FT[0,-1,-2,-3], collect(FT,0:1:20));
roots_flow!(palm.roots, palm.cache_k, palm.cache_p, palm.cache_q, FT(1));
```




## Leaf Hydraulics
The stomatal models often require plant hydraulics either as a correction
    factor (in empirical stomatal models) or as the risk term (in optimal
    stomatal models). To facilitate the calculations, a few specific functions
    are provided.

Function [`xylem_risk`](@ref) returns the risk in xylem hydraulic function
    based on the most downstream end of the xylem. The risk of plant hydraulic
    system is not only on current system, but also potential new growth (plants
    don't want to risk new growth either). Thus, function
    [`xylem_risk`](@ref) evaluates the risk from the xylem pressure
    calculated from current system (with drought history), and then compute the
    risk from the pressure (the severer the srought history, the higher the
    risk):

```@docs
xylem_risk
```

Note that function [`xylem_risk`](@ref) can work on its own without having
    other organ-level components. For example, by changing the `p_ups` of a
    [`LeafHydraulics`](@ref), one can simulate the case of drought without
    caring about other hydraulic systems. Same for function
    [`critical_flow`](@ref) below. However, these functions are only useful for
    sensitivity analysis or when `p_ups` in the [`LeafHydraulics`](@ref) is
    accurate.

Examples

```julia
FT = Float32;
leaf = LeafHydraulics{FT}();
risk = xylem_risk(leaf, FT(0.01));
@show risk;
leaf.p_ups = FT(-1.0);
risk = xylem_risk(leaf, FT(0.01));
@show risk;
```

Function [`critical_flow`](@ref) calculates critical leaf transpiration rate,
    beyond which leaf will desicate. Function [`critical_flow`](@ref) accounts
    for drought legacy effect by design, and the more severe the drought
    history, the lower the `critical_flow`. Again, `ini` in the function is also
    optional, but a good guess will speed up the calculations.

Examples
```julia
FT = Float32;
leaf = LeafHydraulics{FT}();
risk = critical_flow(leaf);
@show risk;
leaf.p_ups = FT(-1.0);
risk = critical_flow(leaf);
@show risk;
```

## Whole-plant Hydraulics
Though [`xylem_risk`](@ref) and [`critical_flow`](@ref) can work on their
    own, the functions only evaluate the risks on leaf level. The more
    realistic case is that when leaf transpiration rate increases, `p_ups` in
    the [`LeafHydraulics`](@ref) gets more negative. Thus, the
    [`xylem_risk`](@ref) and [`critical_flow`](@ref) tends to underestimate
    the risk and overestimate the critical flow rate. To overcome this problem,
    whole-plant level plant hydraulics are provided.

Function [`end_pressure`](@ref) calculates the leaf xylem end pressure for
    a whole-plant struct using these steps:

- calculate the plant base pressure from a given total flow rate
- calculate the trunk end pressure (if present)
- calculate the branch end pressure (if present)
- calculate the leaf end pressure (if present)

Accordingly, there is a function [`critical_flow`](@ref) to calculate the
    critical flow rate for the whole plant. Be aware that Plant-level function
    [`end_pressure`](@ref) and [`critical_flow`](@ref) only applies to the
    case of only one canopy layer (or big-leaf model). As to the case of
    multiple canopy layer, more functions are pending.

```@docs
critical_flow
```

Note that the organ level or whole-plant level conductances are different from
    xylem hydraulic conductance at a given xylem slice. Also, simply computing
    the conductance as the flow rate divided by pressre drop is not accurate
    because of the gravity. In the PlantHydraulics module, the organ and
    whole-plant level conductances are computed by firstly adding up all the
    resistances in each element, and then computing the relative loss of
    conductance. The function available to use is

```@docs
plant_conductances!
```




## Temperature effects
Plant hydraulic properties changes with temperature, due to its impacts on
    surface tension and viscosity. As to surface tension, when temperature
    increases, air-water interface surface tension decreases, meaning that with
    the same curvature, the capillary pressure provided decreases. As a result,
    soil matrix potential becomes less nagative (good for plants!), but the
    conduit resistance to cavitation decreases (bad for plants!). As to
    viscosity, when temperature increases, viscosity decreases, meaning that
    pressure drop decreases at the same flow rate (good for plants!). To
    account for these effects, we provided [`temperature_effects!`](@ref)
    function:

```@docs
temperature_effects!
```

Keep in mind that, leaf critical pressure changes due to surface tension,
    though the impact can be neglected for minor temperature change. However,
    soil water content is still computed using the equivalent matrix potential
    at 25 Celcius because water content is only related to the air-water
    interface curvature.




## Pressure-volume curve

```@docs
AbstractCapacity
PVCurveLinear
PVCurveSegmented
p_from_volume
update_PVF!
```
