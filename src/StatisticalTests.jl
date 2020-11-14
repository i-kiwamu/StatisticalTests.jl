module StatisticalTests

    # packages
    using Printf
    using Distributions, Statistics, StatsBase, GLM
    using Reexport
    @reexport using StatsModels
    # end packages

    # imports
    import Base: show
    import StatsBase: PValue, StatisticalModel, fit, coeftable
    import GLM: LinearModel
    # end imports

    # exports
    export coeftable, fit
    export
        # functions
        t_test
    # end export

    # types
    const NumVector{T<:Number} = AbstractVector{T}
    const IntVector{Int} = AbstractVector{Int}
    const FloatVector{T<:AbstractFloat} = AbstractVector{T}
    abstract type TestModel <: StatisticalModel end
    abstract type TestResult end
    # end types

    # includes
    include("ttest.jl")
    # end includes

end # module
