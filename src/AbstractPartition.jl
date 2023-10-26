abstract type AbstractPartition end

function composition(p::T, q::T) where {T <: AbstractPartition}
    composition_loops(p, q)[1]
end
