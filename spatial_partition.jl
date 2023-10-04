import Base.hash
import Base.==
import Base.copy
using InteractiveUtils
include("partition.jl")

"""
SpatialPartition

Initialize Spatial Partition object

# Arguments
- `upper_points`: upper points of Spatial Partition as array of arrays (length of arrays in arrays are the depth m)
- `lower_points`: lower points of Partition as array of arrays (length of arrays in arrays are the depth m)
- `normal_form`: optional input, transforms to normal form if true, else spatial partition syntax remains
"""
struct SpatialPartition
    upper_points::Array{Int64, 1}
    lower_points::Array{Int64, 1}
    dimension::Int64

    function SpatialPartition(upper_points::AbstractArray, lower_points::AbstractArray, dimension::Int64, normal_form::Bool = true)
        if normal_form
            (upper_points, lower_points) = normal_form_array([Int64.(upper_points), Int64.(lower_points)])
        end
        return new(upper_points, lower_points, dimension)
    end
end

function hash(p::SpatialPartition)

    hash([hash(p.upper_points), hash(p.lower_points), hash(p.dimension)])
    
end

function ==(p::SpatialPartition, q::SpatialPartition)

    p.lower_points == q.lower_points && p.upper_points == q.upper_points && p.dimension == q.dimension

end

function copy(p::SpatialPartition)
    return SpatialPartition(copy(p.upper_points), copy(p.lower_points), copy(p.dimension))
end

"""
helper_new_point_values(x)

This function outputs a semantically identical Partition which has new number values (in O(n)).

# Arguments
- `p`: The input partition which we do not change
- `q`: The input partition that we change according to the values of p

# Returns
- Semantically equal partition to q without point numbers in p
"""
function helper_new_point_values(p::SpatialPartition, q::SpatialPartition)

    p_points::Array{Int} = vcat(copy(p.upper_points), copy(p.lower_points))

    if !(isempty(p_points))
        new_id, index = findmax(p_points)
        new_id += 1
    else
        new_id = 1 
    end
    
    q_points = vcat(copy(q.upper_points), copy(q.lower_points))
    new_ids = Dict()

    for n in q_points
        new_ids[n] = new_id
        new_id += 1
    end

    upper = []
    for (i, n) in enumerate(q.upper_points)
        push!(upper, get(new_ids, n, -1)::Int)
    end

    lower = []
    for (i, n) in enumerate(q.lower_points)
        push!(lower, get(new_ids, n, -1)::Int)
    end

    SpatialPartition(upper, lower, q.dimension, false)
    
end

"""
normal_form(x)

This function outputs a semantically identical Partition which has new number values from 1 to number of blocks of Partitions (in O(n)).
**not needed because we already correct the Partition to the correct form in the constructor**

# Arguments
- `p`: Input partition

# Returns
- `p` with consisten form
"""
function normal_form(p::SpatialPartition)

    new_id::Int = 1
    new_ids::Dict{Int, Int} = Dict()
    p::SpatialPartition = helper_new_point_values(SpatialPartition([length(p.upper_points) + length(p.lower_points)], [], p.dimension), p)

    for (i, n) in enumerate(p.upper_points)
        if !(n in keys(new_ids))
            new_ids[n] = new_id
            p.upper_points[i] = new_id
            new_id += 1
        else
            p.upper_points[i] = get(new_ids, p.upper_points[i], -1)
        end
    end

    for (i, n) in enumerate(p.lower_points)
        if !(n in keys(new_ids))
            new_ids[n] = new_id
            p.lower_points[i] = new_id
            new_id += 1
        else
            p.lower_points[i] = get(new_ids, p.lower_points[i], -1)
        end
    end
    p
end

"""
tensor_product(p, q)

This function applies on p tensor product with q (in O(n)).

# Arguments
- `p`: Input Spatial partition
- `q`: Second input Spatial partition

# Returns
- `p` tensor product `q`
"""
function tensor_product(p::SpatialPartition, q::SpatialPartition)

    @assert p.dimension == q.dimension "p and q have different dimensions in tensor product"

    q_new = helper_new_point_values(p, q)

    SpatialPartition(vcat(p.upper_points, q_new.upper_points), vcat(p.lower_points, q_new.lower_points), p.dimension)
end

"""
involution(p)

This function applies an involution on `p` (in O(n) because normal_form, else O(1)).

# Arguments
- `p`: Input spatial partition

# Returns
- involution of `p`
"""
function involution(p::SpatialPartition)

    SpatialPartition(copy(p.lower_points), copy(p.upper_points), p.dimension)

end

"""
composition(p, q)

This function applies composition between p and q (in O(nlogn)).

# Arguments
- `p`: Input Spatial partition
- `q`: Second input spatial partition

# Returns
- `p` composition `q`
"""
function composition(p::SpatialPartition, q::SpatialPartition)

    @assert p.dimension == q.dimension "p and q have different dimensions in composition"

    comp = composition(Partition(p.upper_points, p.lower_points), Partition(q.upper_points, q.lower_points))

    SpatialPartition(comp.upper_points, comp.lower_points, p.dimension)

end

"""
size(p)

This function outputs the size of the input spatial partition (i.e. the number of points in `p`).

# Arguments
- `p`: Input spatial partition

# Returns
- size of `p`
"""
function size(p::SpatialPartition)
    length(p.lower_points) + length(p.upper_points)
end

function add_partition_to_dict(dict::Dict, p::SpatialPartition)

    add_apbs::Set = get(dict, size(p), -1)
    push!(add_apbs, p)
    dict[size(p)] = add_apbs

    return dict
end

function add_partition_to_composition_dict(array::Array, p::SpatialPartition)

    """add right partition in first dict for top size"""
    add_apbs_top::Set = get(array[1], length(p.upper_points), -1)
    push!(add_apbs_top, p)
    (array[1])[length(p.upper_points)] = add_apbs_top

    """add right partition in first dict for bottom size"""
    add_apbs_bottom::Set = get(array[2], length(p.lower_points), -1)
    push!(add_apbs_bottom, p)
    (array[2])[length(p.lower_points)] = add_apbs_bottom

    return array
end

function do_unary(
    to_unary::Set, 
    all_partitions::Set, 
    stop_whole::Bool, 
    already_u::Set, 
    max_length::Int, 
    all_partitions_by_size::Dict, 
    all_partitions_by_size_top_bottom::Array, 
    trace::Dict)

    stop::Bool = false
    while !stop
        stop = true

        to_unary_copy::Set{SpatialPartition} = copy(to_unary)

        for pp in to_unary_copy
            
            pmod::SpatialPartition = copy(pp)

            a::SpatialPartition = SpatialPartition([], [], 0)

            """continue with involution"""
            a = involution(pmod)

            """add to all_partitions"""
            if !(a in all_partitions)
                trace[a] = tuple([pmod, "i"])
                stop_whole = false
                stop = false
                push!(all_partitions, a)
                push!(to_unary, a)

                """call functions which adds the partition a into the right set in the dict"""
                all_partitions_by_size = add_partition_to_dict(all_partitions_by_size, a)
                all_partitions_by_size_top_bottom = add_partition_to_composition_dict(
                        all_partitions_by_size_top_bottom, a)
            end

            """remember already unary"""
            push!(already_u, pp)
            pop!(to_unary, pp)
        end
    end

    return [stop_whole, all_partitions, already_u, all_partitions_by_size, all_partitions_by_size_top_bottom, trace]
end

function do_tensor_products(
    all_partitions::Set, 
    already_t::Set, 
    to_tens::Set, 
    stop_whole::Bool, 
    max_length::Int, 
    all_partitions_by_size::Dict, 
    all_partitions_by_size_top_bottom::Array, 
    trace::Dict)::Array

    """analogical to all_pyrtitions_by_size in build function for new_tens"""
    new_tens_by_size = Dict()
    for i in 0:max_length
        new_tens_by_size[i] = Set()
    end

    """store all partitions which are new constructed by tensor product"""
    new_tens::Set{SpatialPartition} = Set()

    """store for every i the ii's which are already used, to not use them in this iteration again"""
    without::Dict{SpatialPartition, Set{SpatialPartition}} = Dict()

    """until no more new possibilities tensor"""
    stop::Bool = false
    while !stop
        stop = true

        """if there are new partitions due to tensor and size constraint, remove pair which are already 
        calculated """
        if !isempty(new_tens)
            aa::Set{SpatialPartition} = union(new_tens, all_partitions)
            for i in aa
                """get fitting partitions in advance (improve runtime)"""
                new_tens_temp_tensor::Set{SpatialPartition} = Set()
                for key in keys(new_tens_by_size)
                    if size(i) + Int(key) <= max_length
                        new_tens_temp_tensor = union(new_tens_temp_tensor, get(new_tens_by_size, key, -1)::Set)
                    end
                end
                if i in keys(without)
                    for ii in setdiff(new_tens_temp_tensor, get(without, i, -1))
                        if length(i.upper_points) + length(i.lower_points) + length(ii.upper_points) + length(ii.lower_points) <= max_length && !([i, ii] in already_t)
                            push!(to_tens, [i, ii])
                            push!(already_t, [i, ii])
                        end
                    end
                else
                    for ii in new_tens_temp_tensor
                        if length(i.upper_points) + length(i.lower_points) + length(ii.upper_points) + length(ii.lower_points) <= max_length && !([i, ii] in already_t)
                            push!(to_tens, [i, ii])
                            push!(already_t, [i, ii])
                        end
                    end
                end
            end
        end

        """do the tensor products"""
        al::Set = copy(to_tens)
        for (i, ii) in al
            a::SpatialPartition = tensor_product(i, ii)
            pop!(to_tens, [i, ii])
            if !(a in all_partitions)
                trace[a] = ((i, ii), "t")
                if size(a) == max_length
                    push!(all_partitions, a)

                    """call function which adds the partition a into the right set in the dicts"""
                    all_partitions_by_size = add_partition_to_dict(all_partitions_by_size, a)
                    all_partitions_by_size_top_bottom = add_partition_to_composition_dict(
                        all_partitions_by_size_top_bottom, a)

                    stop_whole = false
                else
                    push!(all_partitions, a)

                    """call function which adds the partition a into the right set in the dicts"""
                    all_partitions_by_size = add_partition_to_dict(all_partitions_by_size, a)
                    all_partitions_by_size_top_bottom = add_partition_to_composition_dict(
                        all_partitions_by_size_top_bottom, a)
                    new_tens_by_size = add_partition_to_dict(new_tens_by_size, a)

                    stop_whole = false
                    push!(new_tens, a)
                    stop = false
                end
            else
                """remove not fitting candidates for further iterations"""
                if !(i in keys(without))
                    without[i] = Set([ii])
                else
                    push!(get(without, i, -1), ii)
                end
            end
        end
    end
    return [all_partitions, already_t, stop_whole, all_partitions_by_size, all_partitions_by_size_top_bottom, trace]
end

function do_composition(
    all_partitions::Set, 
    already_c::Set, 
    stop_whole::Bool, 
    max_length::Int, 
    to_comp::Set, 
    all_partitions_by_size::Dict, 
    all_partitions_by_size_top_bottom::Array, 
    trace::Dict)::Array

    """add newfound partitions due comp"""
    new_comp::Set{SpatialPartition} = Set()

    """new_comp stored in tuple with a dict for top and bottom size (analogical to the technique in build function)"""
    new_comp_by_size_top_bottom = [Dict(), Dict()]
    for i in 0:max_length
        (new_comp_by_size_top_bottom[1])[i] = Set()
        (new_comp_by_size_top_bottom[2])[i] = Set()
    end

    """store for every i the ii's which are already used, to not use them in this iteration again"""
    without::Dict{SpatialPartition, Set{SpatialPartition}} = Dict()

    """until no more new possibilities compose"""
    stop::Bool = false
    while !stop
        stop = true
        """if there are new partitions due to composition, remove pair which are already calculated"""
        if !isempty(new_comp)
            aa = union(new_comp, all_partitions)
            for i in aa
                """get fitting partitions in advance (improve runtime)"""
                new_comp_temp_comp = Set()
                if length(i.upper_points) <= max_length
                    new_comp_temp_comp = get(new_comp_by_size_top_bottom[2], length(i.upper_points), -1)
                end
                if i in keys(without)
                    for ii in setdiff(new_comp_temp_comp, get(without, i, -1))
                        if length(i.upper_points) == length(ii.lower_points) && length(i.upper_points) != 0 && length(i.upper_points) != max_length && length(i.lower_points) + length(ii.upper_points) <= max_length
                            push!(to_comp, [i, ii])
                            push!(already_c, [i, ii])
                        end
                        if length(ii.upper_points) == length(i.lower_points) && length(ii.upper_points) != 0 && length(ii.upper_points) != max_length && length(ii.lower_points) + length(i.upper_points) <= max_length
                            push!(to_comp, [ii, i])
                            push!(already_c, [ii, i])
                        end
                    end
                else
                    for ii in new_comp_temp_comp
                        if length(i.upper_points) == length(ii.lower_points) && length(i.upper_points) != 0 && length(i.upper_points) != max_length && length(i.lower_points) + length(ii.upper_points) <= max_length
                            push!(to_comp, [i, ii])
                            push!(already_c, [i, ii])
                        end
                        if length(ii.upper_points) == length(i.lower_points) && length(ii.upper_points) != 0 && length(ii.upper_points) != max_length && length(ii.lower_points) + length(i.upper_points) <= max_length
                            push!(to_comp, [ii, i])
                            push!(already_c, [ii, i])
                        end
                    end
                end
            end
        end

        """do the compositions"""
        al = copy(to_comp)

        for (i, ii) in al
            a = composition(i, ii)
            if !(a in all_partitions)
                trace[a] = ((i, ii), "c")
                push!(all_partitions, a)

                """call function which adds the partition a into the right set in the dicts"""
                all_partitions_by_size = add_partition_to_dict(all_partitions_by_size, a)
                all_partitions_by_size_top_bottom = add_partition_to_composition_dict(
                    all_partitions_by_size_top_bottom, a)
                new_comp_by_size_top_bottom = add_partition_to_composition_dict(new_comp_by_size_top_bottom, a)

                stop_whole = false
                push!(new_comp, a)
                stop = false
            else
                """remove not fitting candidates for further iterations"""
                pop!(to_comp, [i, ii])
                if !(i in keys(without))
                    without[i] = Set([ii])
                else
                    push!(get(without, i, -1), ii)
                end
            end
        end
    end

    return [all_partitions, already_c, stop_whole, all_partitions_by_size, all_partitions_by_size_top_bottom, trace]
end


"""
construct_category(p::Array, n::Int, tracing::Bool = false, max_artifical::Int = 0)

This function outputs list of all spatial partitions size n constructed from partitions in p (without using spatial partitions of size
greater than max(max(n, max(p)), max_artifical))

# Arguments
- `p`: list of spatial partitions
- `n`: size of spatial partitions in constructing category
- `tracing`: optinal input: activate tracing and get the output (category, trace)
- `max_artifical`: optional input: allow spatial partitions to grow greater while construction process

# Returns
- list of all partitions size n constructed from partitions in p
"""
function construct_category(p::Array, n::Int, tracing::Bool = false, max_artificial::Int = 0)

    """store all candidates found + base partitions to m = 3"""
    all_partitions = Set()
    dim = 0
    if !isempty(p)
        dim = p[1].dimension
        for i in p
            @assert i.dimension == dim "all input partitions need same dimension"
        end
    end
    

    """add spatial base partitions to our set"""
    push!(all_partitions, SpatialPartition([i for i in 1:dim], [i for i in 1:dim], dim))
    push!(all_partitions, SpatialPartition([], vcat([i for i in 1:dim], [i for i in 1:dim]), dim))

    """all candidates stored in dict from size to partition"""
    all_partitions_by_size::Dict{Int, Set{SpatialPartition}} = Dict()

    """all candidates stored in tuple with a dict for top and bottom size"""
    all_partitions_by_size_top_bottom::Array = [Dict(), Dict()]

    """store partitions already unary"""
    already_u::Set{SpatialPartition} = Set()

    """store partitions already tensor product"""
    already_t::Set{Array{SpatialPartition}} = Set()

    """store partitions already composition"""
    already_c::Set{Array{SpatialPartition}} = Set()

    """end output: All partitions found of size n """
    all_partitions_of_size_n::Set{SpatialPartition} = Set()

    """all candidates for unary operations"""
    to_unary::Set{SpatialPartition} = Set(copy(p))

    """trace for tracing"""
    trace::Dict = Dict()

    """compare allowed expansion size with max(n, max_length)"""
    max_length::Int = n

    """get max length of a partition"""
    for i in vcat(p, [SpatialPartition([i for i in 1:dim], [i for i in 1:dim], dim)])
        if size(i) > max_length
            max_length = size(i)
        end
    end

    if max_artificial > 0
        max_length = max(max_length, max_artificial)
    end

    """define for all i <= size an empty set in which we fill the corresponding partition of size i (for tensor)"""
    for i in 0:max_length
        all_partitions_by_size[i] = Set()
    end

    """define for all bottom and top size an empty set in which we fill the corresponding partition"""
    for i in 0:max_length
        (all_partitions_by_size_top_bottom[1])[i] = Set()
        (all_partitions_by_size_top_bottom[2])[i] = Set()
    end

    """add all partitions in p to all_partitions_by_size and all_partitions_by_size_top_bottom"""
    tuple_list_all_partitions::Array = []
    for i in all_partitions
        push!(tuple_list_all_partitions, i)
    end
    for i in vcat(p, tuple_list_all_partitions)
        all_partitions_by_size = add_partition_to_dict(all_partitions_by_size, i)
        all_partitions_by_size_top_bottom = add_partition_to_composition_dict(all_partitions_by_size_top_bottom, i)
    end

    """while new were found apply on them unary tensor and composition"""
    stop_whole::Bool = false
    while !stop_whole
        stop_whole = true
        """add new found partitions in the unary operation candidate list"""
        for i in all_partitions
            if !(i in already_u)
                push!(to_unary, i)
            end
        end

        """fist phase: all possible combinations of unary operations"""
        (stop_whole, all_partitions, already_u, all_partitions_by_size, all_partitions_by_size_top_bottom, trace) = do_unary(to_unary, all_partitions, stop_whole, already_u, max_length, all_partitions_by_size, all_partitions_by_size_top_bottom, trace)

        """store pairs that are candidates to get tensor product"""
        to_tens::Set{Array} = Set()

        """get all pairs to tensor"""
        for i in all_partitions
            """get fitting partitions in advance (improve runtime)"""
            all_partitions_temp_tensor = Set()
            for key in keys(all_partitions_by_size)
                if size(i) + Int(key) <= max_length
                    all_partitions_temp_tensor = union(all_partitions_temp_tensor, get(all_partitions_by_size, key, -1))
                end
            end
            for ii in all_partitions_temp_tensor
                if !((i, ii) in already_t)
                    if length(i.upper_points) + length(i.lower_points) + length(ii.upper_points) + length(ii.lower_points) <= max_length
                        push!(to_tens, [i, ii])
                        push!(already_t, [i, ii])
                    end
                end
            end
        end
        
        """second phase: all possible tensor product operations which aren't redundant (don't do tensor products 
        twice) """
        (all_partitions, already_t, stop_whole, all_partitions_by_size, all_partitions_by_size_top_bottom, trace) = do_tensor_products(all_partitions, already_t, to_tens, stop_whole, max_length, all_partitions_by_size, all_partitions_by_size_top_bottom, trace)

        """add new variations by tensor product or composition with all others"""
        to_comp::Set{Array} = Set()

        """get all pairs to compose"""
        for i in all_partitions
            """get in advance the right second candidate (regarding format)"""
            all_partitions_temp_comp = get(all_partitions_by_size_top_bottom[2], length(i.upper_points), -1)
            for ii in all_partitions_temp_comp
                if !((i, ii) in already_c)
                    if length(i.upper_points) == length(ii.lower_points) && length(i.upper_points) != 0 && length(i.upper_points) != max_length && length(i.lower_points) + length(ii.upper_points) <= max_length
                        push!(to_comp, [i, ii])
                        push!(already_c, [i, ii])
                    end
                end
            end
        end

        """third phase: all possible compositions which aren't redundant (don't do tensor products twice)"""
        (all_partitions, already_c, stop_whole, all_partitions_by_size, all_partitions_by_size_top_bottom, trace) = do_composition(all_partitions, already_c, stop_whole, max_length, to_comp, all_partitions_by_size, all_partitions_by_size_top_bottom, trace)
    end

    """remove all partitions without size n"""
    for i in all_partitions
        if size(i) == n
            if !(i in all_partitions_of_size_n)
                push!(all_partitions_of_size_n, i)
            end
        end
    end

    """format every tuple to partition and return"""

    partitions::Array = []
    for i in all_partitions_of_size_n
        push!(partitions, i)
    end

    if tracing
        return partitions, trace
    end

    return partitions
end

a = SpatialPartition([], [1, 1], 2)
b = SpatialPartition([1, 2, 1, 3], [4, 2, 4, 3], 2)

c = SpatialPartition([1, 2, 3, 2], [1, 4, 3, 4], 2)
dd = SpatialPartition([1, 2, 3, 4], [3, 2, 1, 4], 2)

(p, t) = construct_category([a, c, dd], 8, true, 10)

get_trace(t, b)
println("..............................")
