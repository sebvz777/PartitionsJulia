import Base.hash
import Base.==
import Base.copy
using InteractiveUtils


"""
helper_new_point_values_array(x)

This function outputs a semantically identical Partition in form of an array which has new number values (in O(n)).

# Arguments
- `p`: The input partition as array which we do not change
- `q`: The input partition as array that we change according to the values of p

# Returns
- Semantically equal partition as array to q without point numbers in p
"""
function helper_new_point_values_array(p::Array, q::Array)

    p_points = vcat(copy(p[1]), copy(p[2]))

    if !(isempty(p_points))
        new_id, index = findmax(p_points)
        new_id += 1
    else
        new_id = 1 
    end

    q_points = vcat(q[1], q[2])
    new_ids = Dict()

    for n in q_points
        new_ids[n] = new_id
        new_id += 1
    end

    upper = []
    for (i, n) in enumerate(q[1])
        push!(upper, get(new_ids, n, -1))
    end

    lower = []
    for (i, n) in enumerate(q[2])
        push!(lower, get(new_ids, n, -1))
    end

    [upper, lower]
    
end

"""
normal_form_array(x)

This function outputs a semantically identical Partition of array form which has new number values from 1 to number of blocks of Partitions (in O(n)).

# Arguments
- `p`: Input partition as array

# Returns
- `p` with consisten form
"""
function normal_form_array(p::Array)

    new_id = 1
    new_ids = Dict()
    p = helper_new_point_values_array([[length(p[1]) + length(p[2])], []], p)

    for (i, n) in enumerate(p[1])
        if !(n in keys(new_ids))
            new_ids[n] = new_id
            p[1][i] = new_id
            new_id += 1
        else
            p[1][i] = get(new_ids, p[1][i], -1)
        end
    end

    for (i, n) in enumerate(p[2])
        if !(n in keys(new_ids))
            new_ids[n] = new_id
            p[2][i] = new_id
            new_id += 1
        else
            p[2][i] = get(new_ids, p[2][i], -1)
        end
    end
    p
end

"""
Partition

Initialize Partition object

# Arguments
- `upper_points`: upper points of Partition as array
- `lower_points`: lower points of Partition as array
- `normal_form`: optional input, transforms to normal form if true, else partition syntax remains
"""
struct Partition
    upper_points::Array{Int64, 1}
    lower_points::Array{Int64, 1}

    function Partition(upper_points::AbstractArray, lower_points::AbstractArray, normal_form::Bool = true)
        if normal_form
            (upper_points, lower_points) = normal_form_array([Int64.(upper_points), Int64.(lower_points)])
        end
        return new(upper_points, lower_points)
    end
end


function hash(p::Partition)

    hash([hash(p.upper_points), hash(p.lower_points)])
    
end

function ==(p::Partition, q::Partition)

    p.lower_points == q.lower_points && p.upper_points == q.upper_points

end

function copy(p::Partition)
    return Partition(copy(p.upper_points), copy(p.lower_points))
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
function helper_new_point_values(p::Partition, q::Partition)

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

    Partition(upper, lower, false)
    
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
function normal_form(p::Partition)

    new_id::Int = 1
    new_ids::Dict{Int, Int} = Dict()
    p::Partition = helper_new_point_values(Partition([length(p.upper_points) + length(p.lower_points)], []), p)

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
- `p`: Input partition
- `q`: Second input partition

# Returns
- `p` tensor product `q`

# Examples
```julia-repl
julia> tensor_product(Partition([1, 2], [2, 1]), Partition([1, 1], [1]))
Partition([1, 2, 3, 3], [2, 1, 3])
```
"""
function tensor_product(p::Partition, q::Partition)

    q_new = helper_new_point_values(p, q)

    Partition(vcat(p.upper_points, q_new.upper_points), vcat(p.lower_points, q_new.lower_points))
end

"""
involution(p)

This function applies an involution on `p` (in O(n) because normal_form, else O(1)).

# Arguments
- `p`: Input partition

# Returns
- involution of `p`

# Examples
```julia-repl
julia> involution(Partition([1, 2, 3], [2, 1]))
Partition([1, 2], [2, 1, 3])
```
"""
function involution(p::Partition)

    Partition(copy(p.lower_points), copy(p.upper_points))

end

"""
vertical_reflection(x)

This function applies an vertical reflection on `p` (in O(n) because normal_form, else O(1)).

# Arguments
- `p`: Input partition

# Returns
- vertical reflection of `p`

# Examples
```julia-repl
julia> vertical_reflection(Partition([1, 2, 3], [2, 1]))
Partition([1, 2, 3], [3, 2])
```
"""
function vertical_reflection(p::Partition)

    Partition(reverse(p.upper_points), reverse(p.lower_points))

end

"""
rotation(x)

This function applies a rotation on `p` (in O(n) because normal_form, else O(1)). 

# Arguments
- `p`: Input partition
- `lr`: lr whether left (true) or right (false)
- `tb`: tb whether top (true) or bottom (false) rotation

# Returns
- rotation of `p`

# Examples
```julia-repl
julia> rotation(Partition([1, 2, 3], [2, 1]), true, true)
Partition([1, 2], [3, 1, 3])
```
"""
function rotation(p::Partition, lr::Bool, tb::Bool)

    if tb
        @assert !isempty(p.upper_points) ["Got no partition reaching top"]
    end
    if !tb
        @assert !isempty(p.lower_points) ["Got no partition reaching bottom"]
    end

    ret::Array = [copy(p.upper_points), copy(p.lower_points)]

        if lr
            if tb
                a::Int = ret[1][1]
                splice!(ret[1], 1)
                pushfirst!(ret[2], a)
            else
                a = ret[2][1]
                splice!(ret[2], 1)
                pushfirst!(ret[1], a)
            end
        else
            if tb
                a = ret[1][end]
                pop!(ret[1])
                push!(ret[2], a)
            else
                a = ret[2][end]
                pop!(ret[2])
                push!(ret[1], a)
            end
        end
    Partition(ret[1], ret[2])
end

"""
composition(p, q)

This function applies composition between p and q (in O(nlogn)).

# Arguments
- `p`: Input partition
- `q`: Second input partition
- `loop`: optional input: by default false

# Returns
- `p` composition `q` if loop == false, else [`p` composition `q`, number of loops]

# Examples
```julia-repl
julia> composition(Partition([1, 2], [2, 1]), Partition([1], [1, 1]))
Partition([1], [1, 1])
```
"""
function composition(p::Partition, q::Partition, loop::Bool = false)

    @assert length(p.upper_points) == length(q.lower_points) ["format not fitting"]

    """Work with copies to not change the input partitions"""
    p_copy::Partition = copy(p)

    """new_ids dicts store the new Value we need to assign to the partition in order to connect new segments"""
    q_copy_new_ids::Partition = helper_new_point_values(p_copy, q)
    new_ids::Dict{Int, Int} = Dict()

    """fitting the second partition-values to the first and changing if connection"""
    for (i,n) in enumerate(q_copy_new_ids.lower_points)
        if !(n in keys(new_ids))
            new_ids[n] = p.upper_points[i]
        else
            if p.upper_points[i] in keys(new_ids) && get(new_ids, n, -1) in keys(new_ids)
                """Do path compression if we have the case that we need to merge two tree's together and
                the nodes we operate on are not a root or a leaf"""
                for ii::Int in [n]
                    path::Array = [ii]
                    already_in = Set()
                    push!(already_in, get(new_ids, ii, -1))
                    z::Int = get(new_ids, ii, -1)
                    while z::Int in keys(new_ids)
                        push!(path, z)
                        push!(already_in, z)
                        z = get(new_ids, z, -1)
                        if z in already_in
                            break
                        end
                    end
                    push!(path, z)
                    for nn::Int in path[1:end-1]
                        new_ids[nn] = path[end]
                    end
                end
                new_ids[get(new_ids, n, -1)] = get(new_ids, p.upper_points[i], -1)
            else
                if !(get(new_ids, n, -1) in keys(new_ids))
                    new_ids[get(new_ids, n, -1)] = p.upper_points[i]
                else
                    new_ids[p.upper_points[i]] = get(new_ids, n, 1)
                end
            end
        end
    end
    
    """final path compression"""
    for ii in keys(new_ids)
        path = [ii]
        already_in = Set()
        push!(already_in, get(new_ids, ii, -1))
        z = get(new_ids, ii, -1)
        while z in keys(new_ids)
            push!(path, z)
            push!(already_in, z)
            z = get(new_ids, z, -1)
            if z in already_in
                break
            end
        end
        push!(path, z)
        for nn in path[1:end-1]
            new_ids[nn] = path[end]
        end
    end
    
    """giving the top part new values"""
    for (i,n) in enumerate(q_copy_new_ids.upper_points)
        if n in keys(new_ids)
            q_copy_new_ids.upper_points[i] = get(new_ids, n, -1)
        end
    end

    """giving the top part new values"""
    for (i,n) in enumerate(p_copy.lower_points)
        if n in keys(new_ids)
            p_copy.lower_points[i] = get(new_ids, n, -1)
        end
    end

    """removing the middle by just changing the top of our partition to the adjusted top of the second partition"""
    ret = Partition(q_copy_new_ids.upper_points, p_copy.lower_points)

    """calculating removed related components (loop)"""
        if loop
            related_comp = Set()
            return_partition_as_set = Set(vcat(q_copy_new_ids.upper_points, p_copy.lower_points))

            """calculate new ids for middle nodes, which are under normal circumstances omitted"""
            for (i, n) in enumerate(q_copy_new_ids.lower_points)
                if n in keys(new_ids)
                    q_copy_new_ids.lower_points[i] = get(new_ids, n, -1)
                end
            end

            for (i, n) in enumerate(p_copy.upper_points)
                if n in keys(new_ids)
                    p.upper_points[i] = get(new_ids, n, -1)
                end
            end

            """if there is a ID in the middle part which is not in result partition set we know, that this is a loop"""
            for co in vcat(q_copy_new_ids.lower_points, p_copy.upper_points)
                if !(co in return_partition_as_set)
                    push!(related_comp, co)
                end
            end

            return [ret, length(related_comp)]
        end

    ret
end

"""
size(p)

This function outputs the size of the input partition (i.e. the number of points in `p`).

# Arguments
- `p`: Input partition

# Returns
- size of `p`
"""
function size(p::Partition)
    length(p.lower_points) + length(p.upper_points)
end

function add_partition_to_dict(dict::Dict, p::Partition)::Dict

    add_apbs::Set = get(dict, size(p), -1)
    push!(add_apbs, p)
    dict[size(p)] = add_apbs

    return dict
end

function add_partition_to_composition_dict(array::Array, p::Partition)::Array

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
    trace::Dict)::Array

    stop::Bool = false
    while !stop
        stop = true

        to_unary_copy::Set{Partition} = copy(to_unary)

        for pp in to_unary_copy
            
            pmod::Partition = copy(pp)

            a::Partition = Partition([], [])
            """start with rotation"""
            if !isempty(pmod.upper_points)
                a = rotation(pmod, true, true)
            elseif length(pmod.lower_points) > 0
                a = rotation(pmod, false, false)
            end

            """add to all_partitions"""
            if !(a in all_partitions)
                trace[a] = tuple([pmod, "r"])
                stop_whole = false
                stop = false
                push!(all_partitions, a)
                push!(to_unary, a)

                """call functions which adds the partition a into the right set in the dict"""
                all_partitions_by_size = add_partition_to_dict(all_partitions_by_size, a)
                all_partitions_by_size_top_bottom = add_partition_to_composition_dict(
                        all_partitions_by_size_top_bottom, a)
            end

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

            """end with involution y-axis"""
            a = vertical_reflection(pmod)

            """add to all_partitions"""
            if !(a in all_partitions)
                trace[a] = tuple([pmod, "vr"])
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
    new_tens::Set{Partition} = Set()

    """store for every i the ii's which are already used, to not use them in this iteration again"""
    without::Dict{Partition, Set{Partition}} = Dict()

    """until no more new possibilities tensor"""
    stop::Bool = false
    while !stop
        stop = true

        """if there are new partitions due to tensor and size constraint, remove pair which are already 
        calculated """
        if !isempty(new_tens)
            aa::Set{Partition} = union(new_tens, all_partitions)
            for i in aa
                """get fitting partitions in advance (improve runtime)"""
                new_tens_temp_tensor::Set{Partition} = Set()
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
            a::Partition = tensor_product(i, ii)
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
    new_comp::Set{Partition} = Set()

    """new_comp stored in tuple with a dict for top and bottom size (analogical to the technique in build function)"""
    new_comp_by_size_top_bottom = [Dict(), Dict()]
    for i in 0:max_length
        (new_comp_by_size_top_bottom[1])[i] = Set()
        (new_comp_by_size_top_bottom[2])[i] = Set()
    end

    """store for every i the ii's which are already used, to not use them in this iteration again"""
    without::Dict{Partition, Set{Partition}} = Dict()

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

This function outputs list of all partitions size n constructed from partitions in p (without using partitions of size
greater than max(max(n, max(p)), max_artifical))

# Arguments
- `p`: list of partitions
- `n`: size of partitions in constructing category
- `tracing`: optinal input: activate tracing and get the output (category, trace)
- `max_artifical`: optional input: allow partitions to grow greater while construction process

# Returns
- list of all partitions size n constructed from partitions in p

# Examples
```julia-repl
julia> length(construct_category([Partition([1, 2], [2, 1])], 6))
105
julia> length(construct_category([Partition([1, 2], [2, 1])], 6, true))
[<Partition([1, 2], [2, 1])> ∩ P(6), Dict{Partition -> Tuple}]
```
"""
function construct_category(p::Array, n::Int, tracing::Bool = false, max_artificial::Int = 0)

    """store all candidates found"""
    all_partitions::Set{Partition} = Set([Partition([1, 1], []), Partition([1], [1])])

    """all candidates stored in dict from size to partition"""
    all_partitions_by_size::Dict{Int, Set{Partition}} = Dict()

    """all candidates stored in tuple with a dict for top and bottom size"""
    all_partitions_by_size_top_bottom::Array = [Dict(), Dict()]

    """store partitions already unary"""
    already_u::Set{Partition} = Set()

    """store partitions already tensor product"""
    already_t::Set{Array{Partition}} = Set()

    """store partitions already composition"""
    already_c::Set{Array{Partition}} = Set()

    """end output: All partitions found of size n """
    all_partitions_of_size_n::Set{Partition} = Set()

    """all candidates for unary operations"""
    to_unary::Set{Partition} = Set(copy(p))

    """trace for tracing"""
    trace::Dict = Dict()

    """compare allowed expansion size with max(n, max_length)"""
    max_length::Int = n

    """get max length of a partition"""
    for i in vcat(p, [Partition([1], [1])])
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
                if !((i, ii) in already_t) && (length(i.upper_points) + length(i.lower_points) + length(ii.upper_points) + length(ii.lower_points) <= max_length)
                    push!(to_tens, [i, ii])
                    push!(already_t, [i, ii])
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
                if !((i, ii) in already_c) && (length(i.upper_points) == length(ii.lower_points) && length(i.upper_points) != 0 && length(i.upper_points) != max_length && length(i.lower_points) + length(ii.upper_points) <= max_length)
                    push!(to_comp, [i, ii])
                    push!(already_c, [i, ii])
                end
            end
        end

        """third phase: all possible compositions which aren't redundant (don't do tensor products twice)"""
        (all_partitions, already_c, stop_whole, all_partitions_by_size, all_partitions_by_size_top_bottom, trace) = do_composition(all_partitions, already_c, stop_whole, max_length, to_comp, all_partitions_by_size, all_partitions_by_size_top_bottom, trace)
    end

    """remove all partitions without size n"""
    for i in all_partitions
        if size(i) == n && !(i in all_partitions_of_size_n)
            push!(all_partitions_of_size_n, i)
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

"""
get_trace(trace::Dict, start::Partition)

This function prints out the trace of the partition `start` constructed with `construct_category`
(via breath first search)

"""
function get_trace(trace::Dict, start)
    """track the trace with breath first search"""

    if !(start in keys(trace))
        print("(spatial) Partition $(start) not found in trace")
    end

    track = [start]
    for p in track
        if p in keys(trace)
            println(p, " : ", get(trace, p, -1))
            for i in get(trace, p, -1)[1]
                if !(typeof(i) <: AbstractString) && !(i in track)
                    push!(track, i)
                end
            end
        end
    end
end

"""
check_pair(p::Partition)

This function checks whether `p` is a partition only including blocks of size two (in O(n) average).

# Arguments
- `p`: Input partition

# Returns
- true if `p` is pair partition else false

# Examples
```julia-repl
julia> check_pair(Partition([1, 2, 2, 1, 3], [3]))
true
```
"""
function check_pair(p::Partition)

    """Dictionary from block to size of block"""
    block_to_size = Dict()

    """Initialize dictionary"""
    for i in vcat(p.upper_points, p.lower_points)
        if !(i in keys(block_to_size))
            block_to_size[i] = 1
        else
            block_to_size[i] = get(block_to_size, i, -1) + 1
            if get(block_to_size, i, -1) > 2
                return false
            end
        end
    end
    true
end

"""
check_pair(p::Partition)

This function checks whether `p` is a balanced partition (in O(n) average).

# Arguments
- `p`: Input partition

# Returns
- true if `p` is balanced partition else false

# Examples
```julia-repl
julia> check_balanced(Partition([1, 2, 3], [3, 2, 1]))
true
```
"""
function check_balanced(p::Partition)

    p_array = vcat(p.upper_points, p.lower_points)

    """Dictionary from block to sum of -1 (repr odd indices) and 1 (repr even indices)"""
    block_to_size = Dict()

    """prefill dict with zeros"""
    for i in 1:findmax(p_array)[1]
        block_to_size[i] = 0
    end

    """Initialize dictionary"""
    for (i, n) in enumerate(p_array)
        if i % 2 == 1
            block_to_size[n] = get(block_to_size, n, -1) - 1
        else
            block_to_size[n] = get(block_to_size, n, -1) + 1
        end
    end

    for i in values(block_to_size)
        if i != 0
            return false
        end
    end
    true
end

"""
check_nc(p::Partition)

This function checks whether `p` is a non-crossing partition (in O(n) average).

# Arguments
- `p`: Input partition

# Returns
- true if `p` is non-crossing partition else false

# Examples
```julia-repl
julia> check_nc(Partition([1, 2, 2, 3, 1, 4], [4, 3]))
false
```
"""
function check_nc(p::Partition)

    """transform partition to only upper points"""
    p_array = vcat(p.upper_points, reverse(p.lower_points))

    """Dictionary from block to size of block"""
    block_to_size = Dict()

    """Initialize dictionary"""
    for i in p_array
        if !(i in keys(block_to_size))
            block_to_size[i] = 1
        else
            block_to_size[i] = get(block_to_size, i, -1) + 1
        end
    end

    """blocks we have already seen in the iteration process"""
    already_seen = Set()
    last_incompleted = []
    incompleted = Set()

    for (i, n) in enumerate(p_array)
        if n in incompleted && (isempty(last_incompleted) ? -1 : last_incompleted[end]) != n
            return false
        end
        if !(n in already_seen)
            push!(already_seen, n)
            block_to_size[n] = get(block_to_size, n, -1) - 1
            if get(block_to_size, n, -1) > 0 && (isempty(last_incompleted) ? -1 : last_incompleted[end]) != n
                push!(last_incompleted, n)
                push!(incompleted, n)
            end
        else
            if p_array[i-1] != n && get(block_to_size, p_array[i-1], -1) != 0
                return false
            else
                block_to_size[n] = get(block_to_size, n, -1) - 1
                if get(block_to_size, n, -1) == 0 && !isempty(last_incompleted)
                    if last_incompleted[end] == n
                        pop!(last_incompleted)
                    end
                    delete!(incompleted, n)
                end
                if get(block_to_size, n, -1) > 0 && (isempty(last_incompleted) ? -1 : last_incompleted[end]) != n
                    push!(last_incompleted, n)
                    push!(incompleted, n)
                end
            end
        end
    end
    return true
end
