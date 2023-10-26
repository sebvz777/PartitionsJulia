import Base.hash
import Base.==
import Base.copy
using InteractiveUtils

"""
SpatialPartition

Initialize Spatial Partition object

# Arguments
- `upper_points`: upper points of Spatial Partition as array
- `lower_points`: lower points of Spatial Partition as array 
- `dimension`: dimension of the spatial Partition
- `normal_form`: optional input, transforms to normal form if true, else spatial partition syntax remains
"""
struct SpatialPartition <: AbstractPartition
    partition::Partition
    dimension::Int64
end

function spatial_partition(partition::Partition, dim::Int)

    return SpatialPartition(partition, dim)

end

function hash(p::SpatialPartition, h::UInt)

    hash(p.partition, hash(p.dimension, h))
    
end

function ==(p::SpatialPartition, q::SpatialPartition)

    p.partition == q.partition && p.dimension == q.dimension

end

function copy(p::SpatialPartition)
    return SpatialPartition(copy(p.partition), copy(p.dimension))
end

"""
tensor_product(p::SpatialPartition, q::SpatialPartition)

This function applies on p tensor product with q (in O(n)).

# Arguments
- `p`: Input Spatial partition
- `q`: Second input Spatial partition

# Returns
- `p` tensor product `q`
"""
function tensor_product(p::SpatialPartition, q::SpatialPartition)

    @assert p.dimension == q.dimension "p and q have different dimensions in tensor product"

    SpatialPartition(tensor_product(p.partition, q.partition), p.dimension)
end

"""
involution(p::SpatialPartition)

This function applies an involution on `p` (in O(n) because normal_form, else O(1)).

# Arguments
- `p`: Input spatial partition

# Returns
- involution of `p`
"""
function involution(p::SpatialPartition)

    SpatialPartition(involution(p.partition), p.dimension)

end

"""
composition_loops(p::SpatialPartition, q::SpatialPartition)

This function applies composition between p and q (in O(nlogn)).

# Arguments
- `p`: Input partition
- `q`: Second input partition

# Returns
- [`p` composition `q`, number of loops]
"""
function composition_loops(p::SpatialPartition, q::SpatialPartition)

    @assert is_composable(p, q) "p and q have different dimensions in composition"

    comp_loops = composition_loops(p.partition, q.partition)

    (SpatialPartition(comp_loops[1], p.dimension), comp_loops[2])

end