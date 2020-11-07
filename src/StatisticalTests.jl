module StatisticalTests

    # packages
    using Printf
    using Distributions, Statistics, StatsBase
    # using DataFrames
    # end packages

    # imports
    import Base: show
    import StatsBase: PValue
    # end imports

    # exports
    export
        # functions
        t_test
    # end export

    # types
    const NumVector{T<:Number} = AbstractVector{T}
    const IntVector{Int} = AbstractVector{Int}
    const FloatVector{T<:AbstractFloat} = AbstractVector{T}
    abstract type TestResult end
    # end types

    # includes
    include("ttest.jl")
    # end includes

end # module
