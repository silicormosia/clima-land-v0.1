###############################################################################
#
# Calculate the product-limited photosynthetic rates
#
###############################################################################
"""
    product_limited_rate!(
                photo_set::C3ParaSet{FT},
                leaf::Leaf{FT}
    ) where {FT<:AbstractFloat}
    product_limited_rate!(
                photo_set::C4ParaSet{FT},
                leaf::Leaf{FT}
    ) where {FT<:AbstractFloat}
    product_limited_rate!(
                photo_set::C4ParaSet{FT},
                leaf::Leaf{FT},
                envir::AirLayer{FT}
    ) where {FT<:AbstractFloat}

Calculate the product limited photosynthetic rate, given
- `photo_set` [`C3ParaSet`](@ref) or [`C4ParaSet`](@ref) type struct
- `leaf` [`Leaf`](@ref) type struct
- `envir` [`AirLayer`](@ref) type struct
"""
function product_limited_rate!(
            photo_set::C3Cytochrome{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    (; p_i, Γ_star) = leaf;
    (; Eff_1, Eff_2) = photo_set;

    leaf.Ap = leaf.Vcmax / 2;
    leaf.J_P680_p = leaf.Ac * (Eff_1*p_i + Eff_2*Γ_star) / (p_i - Γ_star);
    leaf.J_P700_p = leaf.J_P680_c * leaf.η;

    return nothing
end




function product_limited_rate!(
            photo_set::C3ParaSet{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    leaf.Ap = leaf.Vcmax / 2;

    return nothing
end




function product_limited_rate!(
            photo_set::C4ParaSet{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    leaf.Ap = leaf.Kpep * leaf.Vcmax25 * leaf.p_i;

    return nothing
end




function product_limited_rate!(
            photo_set::C4ParaSet{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT}
) where {FT<:AbstractFloat}
    (; g_lc, Kpep, Rd, Vcmax25) = leaf;
    (; p_a, p_atm) = envir;

    _p_i = (FT(1e6) * g_lc * p_a + p_atm * Rd) / (p_atm * Kpep * Vcmax25 + FT(1e6) * g_lc);
    leaf.Ap = Kpep * Vcmax25 * _p_i;

    return nothing
end
