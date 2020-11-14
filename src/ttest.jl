mutable struct TTestResult<:TestResult
    type::Symbol  # either :one, :paired, :simple, or :welch
    n_obs::IntVector
    means::FloatVector
    sds::FloatVector
    statistic::Float64
    df::Union{Int, Float64}
    pval::Float64
    cohen_d::Union{Float64, Nothing}
end

# function show(io::IO, ttr::TTestResult)
#     println(coeftable(ttr))
# end

function coeftable(ttr::TTestResult)
    if ttr.type == :one
        println("\n        One sample t-test\n")
        println("Data (mean ± sd)")
        @printf("  x: %.2f ± %.2f\n", ttr.means[1], ttr.sds[1])
        @printf("  difference from %.2f\n", ttr.means[2])
        println("\nResult")
        CoefTable(
            [ttr.statistic, ttr.df, ttr.pval],
            ["t", "df", "P-value"],
            [""],
            3, 1
        )
    elseif ttr.type == :paired
        println("\n        Paired t-test\n")
        @printf("Difference (mean ± sd): %.2f ± %.2f\n",
                ttr.means[1], ttr.sds[1])
        println("\nResult")
        CoefTable(
            [ttr.statistic, ttr.df, ttr.pval],
            ["t", "df", "P-value"],
            [""],
            3, 1
        )
    elseif ttr.type == :simple
        println("\n        Simple Two Sample t-test\n")
        println("Data (mean ± sd)")
        @printf("  x1: %.2f ± %.2f\n", ttr.means[1], ttr.sds[1])
        @printf("  x2: %.2f ± %.2f\n", ttr.means[2], ttr.sds[2])
        println("\nResult")
        CoefTable(
            [ttr.statistic, ttr.df, ttr.pval, ttr.cohen_d],
            ["t", "df", "P-value", "Cohen's d"],
            [""],
            3, 1
        )
    elseif ttr.type == :welch
        println("\n        Welch Two Sample t-test\n")
        println("Data (mean ± sd)")
        @printf("  x1: %.2f ± %.2f\n", ttr.means[1], ttr.sds[1])
        @printf("  x2: %.2f ± %.2f\n", ttr.means[2], ttr.sds[2])
        println("\nResult")
        CoefTable(
            [ttr.statistic, ttr.df, ttr.pval, ttr.cohen_d],
            ["t", "df", "P-value", "Cohen's d"],
            [""],
            3, 1
        )
    else
        error("invalid t-test type $ttr.type")
    end
end

"""
    t_test(x, μ0::AbstractFloat=0.0)

One sample t-test

The argument `x` can be numeric `Vector`.

One sample t-test will be performed to compare between the mean of `x1` and `μ0`.
"""
function t_test(x::NumVector; μ0::AbstractFloat=0.0)::TTestResult
    # check types
    if !(x isa FloatVector)
        @warn "x is converted to Float64"
        x = convert(Vector{Float64}, x)
    end
    
    # check length
    n = length(x)
    n > 1 || error("size of x is not enough: actual size is $n")

    μ = mean(x)
    σ = std(x)
    t_stat = (μ - μ0) / (σ / sqrt(n))
    ν = n - 1
    pval = 2.0 * ccdf(TDist(ν), abs(t_stat))

    res = TTestResult(:one, [n], [μ, μ0], [σ], t_stat, ν, pval, nothing)
    return res
end

"""
    t_test(x1, x2, paired::Bool=false, varequal::Bool=false)

Two sample t-test

The arguments `x1` and `x2` can be numeric `Vector`.

Paired t-test will be performed if the keyword argument `paired` is `true`. Simple two sample t-test will be performed if the keyword argument `varequal` is `true`. Otherwise, Welch two sample t-test will be performed.
"""
function t_test(
    x1::NumVector,
    x2::NumVector;
    paired::Bool=false,
    varequal::Bool=false
)::TTestResult
    # check types
    if !(x1 isa FloatVector)
        @warn "x1 is converted to Float64"
        x1 = convert(Vector{Float64}, x1)
    end
    if !(x2 isa FloatVector)
        @warn "x2 is converted to Float64"
        x2 = convert(Vector{Float64}, x2)
    end

    # check length
    n1 = length(x1)
    n2 = length(x2)
    n1 > 1 || error("size of x1 is not enough: actual size is $n1")
    n2 > 1 || error("size of x2 is not enough: actual size is $n2")

    # paired, simple, or Welch
    if paired
        n1 == n2 || error("size of x1 and x2 must be the same: ($n1, $n2)")
        t_test_paired(x1, x2)
    elseif varequal
        t_test_var_equal(x1, x2)
    else
        t_test_welch(x1, x2)
    end
end

function t_test_paired(x1::FloatVector, x2::FloatVector)::TTestResult
    n = length(x1)
    Δ = x1 .- x2
    μ = mean(Δ)
    σ = std(Δ)
    t_stat = μ / (σ / sqrt(n))
    ν = n - 1
    pval = 2.0 * ccdf(TDist(ν), abs(t_stat))

    res = TTestResult(:two, [n], [μ], [σ], t_stat, ν, pval, nothing)
    return res
end

function t_test_var_equal(x1::FloatVector, x2::FloatVector)::TTestResult
    m1 = mean(x1)
    m2 = mean(x2)
    v1 = var(x1)
    v2 = var(x2)
    n1 = length(x1)
    n2 = length(x2)

    Δ = m1 - m2
    ν = n1 + n2 - 2
    σ = sqrt(((n1-1)*v1 + (n2-1)*v2) / ν)
    t_stat = Δ / (σ * sqrt(1/n1 + 1/n2))
    pval = 2.0 * ccdf(TDist(ν), abs(t_stat))
    cohen_d = abs(Δ) / sqrt((v1 + v2) / 2.0)

    res = TTestResult(:simple, [n1, n2], [m1, m2], [sqrt(v1), sqrt(v2)], 
                      t_stat, ν, pval, cohen_d)
    return res
end

function t_test_welch(x1::FloatVector, x2::FloatVector)::TTestResult
    m1 = mean(x1)
    m2 = mean(x2)
    v1 = var(x1)
    v2 = var(x2)
    n1 = length(x1)
    n2 = length(x2)
    vn1 = v1 / n1
    vn2 = v2 / n2
    
    Δ = m1 - m2
    σ = sqrt(vn1 + vn2)
    t_stat = Δ / σ
    ν = (vn1 + vn2)^2 / ((vn1^2/(n1-1)) + (vn2^2/(n2-1)))
    pval = 2.0 * ccdf(TDist(ν), abs(t_stat))
    cohen_d = abs(Δ) / sqrt((v1 + v2) / 2.0)
    
    res = TTestResult(:welch, [n1, n2], [m1, m2], [sqrt(v1), sqrt(v2)],
                      t_stat, ν, pval, cohen_d)
    return res
end

"""
    TTestModel

A `TestModel` type

# Members

- `X`: model matrix
- `y`: vector to compare its means
- `type`: either :one, :paired, :simple, or :welch
- `unique_levels`: `String` vector of the unique levels (only 1 or 2 elements are allowed)
"""
struct TTestModel{T<:Real} <: TestModel
    x1::AbstractVector{T}
    x2::AbstractVector{T}
    type::Symbol
    unique_levels::Vector{String}
    function TTestModel(
        X::Matrix{T},
        y::AbstractVector{T},
        type::Symbol,
        unique_levels::Vector{String};
        μ0::T = 0.0
    ) where T
        if type == :one
            new{T}(y, [μ0], :one, unique_levels)
        else
            dummy = X[:,2]
            x1 = y[dummy .== 0]
            x2 = y[dummy .== 1]
            new{T}(x1, x2, type, unique_levels)
        end
    end
end

function fit(obj::TTestModel)
    if obj.type == :one
        t_test(obj.x1, μ0=obj.x2[1])
    else
        t_test(obj.x1, obj.x2, paired=(obj.type == :paired), varequal=(obj.type == :simple))
    end
end

function fit(
        ::Type{TTestModel},
        X::AbstractMatrix{<:Real},
        y::AbstractVector{<:Real},
        type::Symbol,
        unique_levels::Vector{String};
        μ0::Real = 0.0
    )
    ttest_obj = TTestModel(X, y, type, unique_levels, μ0=μ0)
    fit(ttest_obj)
end
function fit(
        ::Type{TTestModel},
        X::AbstractMatrix{<:Real},
        y::AbstractVector{<:Real},
        type::Symbol
    )
    fit(TTestModel, X, y, type, ["x1", "x2"], μ0=0.0)
end

"""
    t_test(X, y, type, unique_levels; μ0=0.0)

An alias for `fit(TTestModel, X, y, unique_levels)`

The arguments `X` and `y` can be a `Matrix` and a `Vector` or a `Formula` and `DataFrame`.

The argument `type` must be either `:one`, `:paired`, `:simple`, or `welch`.

The argument `unique_levels` can be a `Vector{String}`.

The keyword argument `μ0` can be a `Real` to use one sample t-test (default = 0.0). Ignored if the `type` is not `:one`.
"""
function t_test(X, y, paired::Bool=false, varequal::Bool=false)
    if paired
        fit(TTestModel, X, y, :paired, ["x1", "x2"])
    elseif varequal
        fit(TTestModel, X, y, :simple, ["x1", "x2"])
    else
        fit(TTestModel, X, y, :welch, ["x1", "x2"])
    end
end
