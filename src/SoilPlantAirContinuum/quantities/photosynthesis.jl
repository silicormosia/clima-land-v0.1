"""

    A_GROSS(node::SPACMono{FT}) where {FT}

Return the gross photosynthesis rate per layer per ground area, given
- `node` [`SPACMono`](@ref) type struct

"""
A_GROSS(node::SPACMono{FT}) where {FT} = [FT[iPS.Ag' * iPS.LAIx * iPS.LA for iPS in node.plant_ps],
                                          FT[iPS.Ac' * iPS.LAIx * iPS.LA for iPS in node.plant_ps],
                                          FT[iPS.Aj' * iPS.LAIx * iPS.LA for iPS in node.plant_ps],
                                          FT[iPS.Ap' * iPS.LAIx * iPS.LA for iPS in node.plant_ps]] ./ node.ga;


"""

    A_GROSS_RD(node::SPACMono{FT}) where {FT}

Return the gross photosynthesis rate per layer per ground area, given
- `node` [`SPACMono`](@ref) type struct

"""
A_GROSS_RD(node::SPACMono{FT}) where {FT} = [FT[max(0, (iPS.Ag .- iPS.ps.Rd)' * iPS.LAIx * iPS.LA) for iPS in node.plant_ps],
                                             FT[max(0, (iPS.Ac .- iPS.ps.Rd)' * iPS.LAIx * iPS.LA) for iPS in node.plant_ps],
                                             FT[max(0, (iPS.Aj .- iPS.ps.Rd)' * iPS.LAIx * iPS.LA) for iPS in node.plant_ps],
                                             FT[max(0, (iPS.Ap .- iPS.ps.Rd)' * iPS.LAIx * iPS.LA) for iPS in node.plant_ps]] ./ node.ga;


"""
    CNPP(spac::SPACMono{FT}) where {FT<:AbstractFloat}

Return the canopy NPP of the SPAC per ground area
"""
function CNPP(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    cnpp::FT = 0;

    for iPS in spac.plant_ps
        cnpp += numerical∫(iPS.An, iPS.LAIx) * iPS.LA;
    end;

    return cnpp / spac.ga
end;


"""
    GPP(spac::SPACMono{FT}) where {FT<:AbstractFloat}

Return the GPP of the SPAC per ground area
"""
function GPP(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    gpp::FT = 0;

    for iPS in spac.plant_ps
        gpp += numerical∫(iPS.Ag, iPS.LAIx) * iPS.LA;
    end;

    return gpp / spac.ga
end;


"""
    GPP_RD(spac::SPACMono{FT}) where {FT<:AbstractFloat}

Return the max(0, GPP - Rd) of the SPAC per ground area
"""
function GPP_RD(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    gpp::FT = 0;

    for iPS in spac.plant_ps
        gpp += max(0, numerical∫(iPS.An, iPS.LAIx) * iPS.LA);
    end;

    return gpp / spac.ga
end;
