"""
    isna(x)

Indicate whether `x` is `missing` or `NaN`.
"""
isna(x::Any) = isnan(x) ? true : false
isna(::Missing) = true


# """
#     skipnas(itr)

# Return an iterator over the elements in `itr` skipping all non-applicable values (`missing` and `NaN`). The returned object can be indexed using indices of `itr`.
# """
# skipnas(itr) = SkipNAs(itr)

# struct SkipNAs{T}
#     x::T
# end
# IteratorSize(::Type{<:SkipNAs}) = SizeUnknown()
# IteratorEltype(::Type{SkipNAs{T}}) where {T} = IteratorEltype(T)
# eltype(::Type{SkipNAs{T}}) where {T} = nonmissingtype(eltype(T))

# function iterate(itr::SkipNAs, state...)
#     item, state = iterate(itr.x, state...)
#     while isna(item)
#         item, state = iterate(itr.x, state)
#     end
#     item, state
# end

# IndexStyle(::Type{<:SkipNAs{T}}) where {T} = IndexStyle(T)
# eachindex(itr::SkipNAs) =
#     Iterators.filter(i -> !isna(@inbounds(itr.x[i])), eachindex(itr.x))
# keys(itr::SkipNAs) =
#     Iterators.filter(i -> !isna(@inbounds(itr.x[i])), keys(itr.x))
# Base.@propagate_inbounds function getindex(itr::SkipNAs, I...)
#     v = itr.x[I...]
#     isna(v) && throw(MissingException("the value at index $I is missing or NaN"))
#     v
# end


"""
    lchoose(n, j)
Return the log of the number of choose j from n
"""
function lchoose(n, j)
    return loggamma(n+1) - loggamma(j+1) - loggamma(n-j+1)
end
