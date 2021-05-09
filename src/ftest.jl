mutable struct FTestResult <: TestResult
    levels::Union{Vector{String}, Nothing}
    n_obs::IntVector
    sds::FloatVector
    statistic::Float64
    dofs::IntVector
    pval::Float64
end

function show(io::IO, ftr::FTestResult)
    println("\n      F test to compare two variances")
    println("Data (± SD)")
    @printf(io, "  %s: ± %.2f\n",
            ftr.levels[1], ftr.sds[1])
    @printf(io, "  %s: ± %.2f\n",
            ftr.levels[2], ftr.sds[2])
    println(io, "\nResult:\n", coeftable(ftr))
end

function coeftable(ftr::FTestResult)
    CoefTable(
        [ftr.statistic, ftr.dofs[1], ftr.dofs[2], ftr.pval],
        ["F", "DF₁", "DF₂", "P-value"],
        [""],
        4, 1
    )
end

"""
    f_test(x1, x2)

F-test

The arguments `x1` and `x2` can be numeric `Vector`.
"""
function f_test(
    x1::NumVector,
    x2::NumVector;
    levels::Union{Vector{String}, Nothing}=nothing
)::FTestResult
    # use temporal sample names
    if isnothing(levels)
        @info "Using the first and second sample names as 'Sample1' and 'Sample2'"
        levels = ["Sample1", "Sample2"]
    end

    # check types
    if !(x1 isa FloatVector)
        @warn "$levels[1] is converted to Float64"
        x1 = convert(Vector{Float64}, x1)
    end
    if !(x2 isa FloatVector)
        @warn "$levels[2] is converted to Float64"
        x2 = convert(Vector{Float64}, x2)
    end

    # check length
    n1 = length(x1)
    n2 = length(x2)
    n1 > 1 || error("size of $levels[1] is not enough: actual size is $n1")
    n2 > 1 || error("size of $levels[2] is not enough: actual size is $n2")

    # perform f-test
    v1 = var(x1)
    v2 = var(x2)
    ν1 = v1 < v2 ? n2-1 : n1-1
    ν2 = v1 < v2 ? n1-1 : n2-1
    f_stat = v1 < v2 ? v2/v1 : v1/v2
    pval = 2.0 * D.ccdf(FDist(ν1, ν2), f_stat)
    res = FTestResult(levels, [n1, n2], [sqrt(v1), sqrt(v2)],
                      f_stat, [ν1, ν2], pval)
    return res
end

"""
    FTestModel

A `TestModel` type for f-test

# Members

- `X`: model matrix
- `y`: vector to compare its means
- `levels`: `String` vector of the unique levels (only 2 elements are allowed)

"""
struct FTestModel{T<:Real} <: TestModel
    x1::AbstractVector{T}
    x2::AbstractVector{T}
    levels::Vector{String}
    function FTestModel(
        X::Matrix{T},
        y::AbstractVector{T},
        levels::Vector{String};
    ) where T
        dummy = X[:,2]
        x1 = y[dummy .== 0]
        x2 = y[dummy .== 1]
        new{T}(x1, x2, levels)
    end
end

function fit(obj::FTestModel)
    f_test(obj.x1, obj.x2, levels=obj.levels)
end
function fit(
    ::Type{FTestModel},
    X::AbstractMatrix{<:Real},
    y::AbstractVector{<:Real},
    levels::Vector{String}
)
    ftest_obj = FTestModel(X, y, levels)
    fit(ftest_obj)
end

"""
    f_test(X, y)

An alias for `fit(FTestModel, X, y, levels)`

The arguments `X` and `y` can be a `Matrix` and a `Vector` or a `Formula` and `DataFrame`.

The argument `levels` can be a `Vector{String}`.
"""
function f_test(
    X::AbstractMatrix{<:Real},
    y::AbstractVector{<:Real}
)
    @info "Using the first and second sample names as 'Sample1' and 'Sample2'"
    fit(FTestModel, X, y, ["Sample1", "Sample2"])
end
function f_test(
    formula::FormulaTerm{Term,Term},
    df::DataFrame
)
    # obtain levels
    f = apply_schema(formula, schema(formula, df))
    levels = DataAPI.levels(df[!,f.rhs.terms[1].sym])
    if eltype(levels) <: Real
        levels = string.(levels)
    end

    # model matrix
    model_frame = ModelFrame(formula, df)
    X = ModelMatrix(model_frame).m
    y = response(model_frame)

    # perform f-test
    fit(FTestModel, X, y, levels)
end
