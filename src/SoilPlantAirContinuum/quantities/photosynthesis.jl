"""

    ag_per_layer(node::SPACMono{FT}) where {FT}

Return the gross photosynthesis rate per layer, given
- `node` [`SPACMono`](@ref) type struct

"""
ag_per_layer(node::SPACMono{FT}) where {FT} = FT[iPS.Ag' * iPS.LAIx for iPS in node.plant_ps];
