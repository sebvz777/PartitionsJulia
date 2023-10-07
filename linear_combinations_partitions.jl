import Base.+
import Base.*
import Base.hash
import Base.==
import Base.copy
using InteractiveUtils
include("partition.jl")

mutable struct LinearCombPartition

    sum::Array{Array}

    function LinearCombPartition(sum::AbstractArray)
        obj = new(sum)
        obj.sum = simplify_term(obj)
        return obj
    end

    function simplify_term(term::LinearCombPartition)
        """
        Simplify term by removing zero summand and distributivity
        """
            partitions = Dict()
    
            for (i, (i1, i2)) in enumerate(copy(term.sum))
                """removing zero summand"""
                if i1 == 0
                    pop!(term.sum, [i1, i2])
                    continue
                end
                """simplifying distributivity"""
                if !(i2 in keys(partitions))
                    partitions[i2] = i1
                else
                    deleteat!(term.sum, findfirst(x -> x == [i1, i2], term.sum))
                    deleteat!(term.sum, findfirst(x -> x == [get(partitions, i2, -1), i2], term.sum))
                    push!(term.sum, [get(partitions, i2, -1) + i1, i2])
                    partitions[i2] = get(partitions, i2, -1) + i1
                end
            end
            term.sum
        end
end

function simplify_term(term::LinearCombPartition)
    """
    Simplify term by removing zero summand and distributivity
    """
        partitions = Dict()

        for (i, (i1, i2)) in enumerate(copy(term.sum))
            """removing zero summand"""
            if i1 == 0
                pop!(term.sum, [i1, i2])
                continue
            end
            """simplifying distributivity"""
            if !(i2 in keys(partitions))
                partitions[i2] = i1
            else
                deleteat!(term.sum, findfirst(x -> x == [i1, i2], term.sum))
                deleteat!(term.sum, findfirst(x -> x == [get(partitions, i2, -1), i2], term.sum))
                push!(term.sum, [get(partitions, i2, -1) + i1, i2])
                partitions[i2] = get(partitions, i2, -1) + i1
            end
        end
        term
    end

function hash(p::LinearCombPartition)

    hash(p.sum)
    
end

function ==(p::LinearCombPartition, q::LinearCombPartition)

    p = simplify_term(p)
    q = simplify_term(q)

    set_q = Set(q)

    for i in p
        if !(i in set_q)
            return false
        end
    end

    return true
end

function copy(p::LinearCombPartition)
    return LinearCombPartition(copy(p.sum))
end

function +(s1::LinearCombPartition, s2::LinearCombPartition)

    """init dict because of time complexity improvements and simplify to prevent information loss"""
    s1 = simplify_term(s1)
    dict_term = Dict(item[2] => item[1] for item in s1.sum)

    """iterate over the objects lists"""
    for i in s2.sum
        if i[2] in keys(dict_term)
            d_term = get(dict_term, i[2], -1)
            d_term += i[1]
            dict_term[i[2]] = d_term
        else
            dict_term[i[2]] = i[1]
        end
    end

    """transform dict to list and output object"""
    out = [[pair[2], pair[1]] for pair in pairs(dict_term)]

    LinearCombPartition(out)
end

function *(s1, s2)

    out = []
    """iterate over the objects lists and do the product for each pair"""
    for i in s1.sum
        for ii in s2.sum
            """get the composition with the loop level"""
            comp, loop = composition(i[2], ii[2], true)
            """add the simplified form to the output with d^loops"""
            push!(out, [i[1]*ii[1]*(1^loop), comp])
        end
    end
    
    LinearCombPartition(out)
end

println(LinearCombPartition([[3, Partition([1], [1])], [1, Partition([1], [1])], [2, Partition([1], [1])]]))
