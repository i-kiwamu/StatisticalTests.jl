module StatisticalTests

    # packages
    using Printf
    using DataAPI, DataFrames
    using Distributions, Statistics, StatsBase, StatsModels, GLM
    using Reexport
    @reexport using StatsModels
    # end packages

    # imports
    import Base: show
    import StatsBase: PValue, StatisticalModel, fit, coeftable
    import DataAPI: levels
    import DataFrames: DataFrame
    import GLM: LinearModel
    import SpecialFunctions: loggamma
    # end imports

    # exports
    export coeftable, fit
    export
        # functions
        t_test,
        f_test,
        ks_test
    # end export

    # types
    const NumVector{T<:Number} = AbstractVector{T}
    const IntVector{Int} = AbstractVector{Int}
    const FloatVector{T<:AbstractFloat} = AbstractVector{T}
    abstract type TestModel <: StatisticalModel end
    abstract type TestResult end
    # end types

    # includes
    include("utilities.jl")
    include("ttest.jl")
    include("ftest.jl")
    include("kstest.jl")
    # end includes

end # module
