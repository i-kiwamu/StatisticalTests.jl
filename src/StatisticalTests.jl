module StatisticalTests

using Printf
using Distributions, Statistics
# using DataFrames

# exports
export
    # functions
    t_test
# end export

# types
const NumVector{T<:Number} = AbstractVector{T}
const FloatVector{T<:AbstractFloat} = AbstractVector{T}
# end types

# includes
include("ttest.jl")
# end includes

end # module
