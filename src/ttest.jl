mutable struct TTestResult<:TestResult
    type::Symbol  # either :one, :paired, :simple, or :welch
    levels::Union{Vector{String}, Nothing}
    n_obs::IntVector
    means::FloatVector
    sds::FloatVector
    ci::FloatMatrix
    statistic::Float64
    dof::Union{Int, Float64}
    pval::Float64
    cohen_d::Union{Float64, Nothing}
end

function show(io::IO, ttr::TTestResult)
    if ttr.type == :one
        println(io, "\n        One Sample t-test\n")
        @printf(io, "Data (Mean ± SD [95%% CI]): %.2f ± %.2f [%.2f, %.2f]\n",
                ttr.means[1], ttr.sds[1], ttr.ci[1,1], ttr.ci[1,2])
        @printf(io, "  difference from %.2f\n",
                ttr.means[2])
    elseif ttr.type == :paired
        println(io, "\n        Paired t-test\n")
        @printf(io, "Difference (Mean ± SD [95%% CI]): %.2f ± %.2f [%.2f, %.2f]\n",
                ttr.means[1], ttr.sds[1], ttr.ci[1,1], ttr.ci[1,2])
    elseif ttr.type == :simple
        println(io, "\n        Simple Two Sample t-test\n")
        println(io, "Data (Mean ± SD [95%% CI])")
        @printf(io, "  %s: %.2f ± %.2f [%.2f, %.2f]\n",
                ttr.levels[1], ttr.means[1], ttr.sds[1],
                ttr.ci[1,1], ttr.ci[1,2])
        @printf(io, "  %s: %.2f ± %.2f [%.2f, %.2f]\n",
                ttr.levels[2], ttr.means[2], ttr.sds[2],
                ttr.ci[2,1], ttr.ci[2,2])
    elseif ttr.type == :welch
        println(io, "\n        Welch Two Sample t-test\n")
        println(io, "Data (Mean ± SD [95%% CI])")
        @printf(io, "  %s: %.2f ± %.2f [%.2f, %.2f]\n",
                ttr.levels[1], ttr.means[1], ttr.sds[1],
                ttr.ci[1,1], ttr.ci[1,2])
        @printf(io, "  %s: %.2f ± %.2f [%.2f, %.2f]\n",
                ttr.levels[2], ttr.means[2], ttr.sds[2],
                ttr.ci[2,1], ttr.ci[2,2])
    end
    println(io, "\nResult:\n", coeftable(ttr))
end

function coeftable(ttr::TTestResult)
    if (ttr.type == :one) || (ttr.type == :paired)
        CoefTable(
            [ttr.statistic, ttr.dof, ttr.pval],
            ["t", "DF", "P-value"],
            [""],
            3, 1
        )
    elseif (ttr.type == :simple) || (ttr.type == :welch)
        CoefTable(
            [ttr.statistic, ttr.dof, ttr.pval, ttr.cohen_d],
            ["t", "DF", "P-value", "Cohen's d"],
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
    distr_t = TDist(ν)
    pval = 2.0 * D.ccdf(distr_t, abs(t_stat))
    ci = μ .+ quantile.(distr_t, [0.025, 0.975]) .* (σ/sqrt(n))

    res = TTestResult(
        :one, nothing, [n], [μ, μ0], [σ], reshape(ci, (1, 2)),
        t_stat, ν, pval, nothing
    )
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
    levels::Union{Vector{String}, Nothing}=nothing,
    paired::Bool=false,
    varequal::Bool=false
)::TTestResult
    # use temporal sample names
    if !paired & isnothing(levels)
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

    # paired, simple, or Welch
    if paired
        n1 == n2 || error("size of $levels[1] and $levels[2] must be the same: ($n1, $n2)")
        t_test_paired(x1, x2)
    elseif varequal
        t_test_var_equal(x1, x2, levels)
    else
        t_test_welch(x1, x2, levels)
    end
end

function t_test_paired(
    x1::FloatVector,
    x2::FloatVector,
)::TTestResult
    n = length(x1)
    Δ = x1 .- x2
    μ = mean(Δ)
    σ = std(Δ)
    t_stat = μ / (σ / sqrt(n))
    ν = n - 1
    distr_t = TDist(ν)
    pval = 2.0 * D.ccdf(distr_t, abs(t_stat))
    ci = μ .+ quantile.(distr_t, [0.025, 0.975]) .* (σ/sqrt(n))

    res = TTestResult(
        :paired, nothing, [n], [μ], [σ], reshape(ci, (1,2)),
        t_stat, ν, pval, nothing
    )
    return res
end

function t_test_var_equal(
    x1::FloatVector,
    x2::FloatVector,
    levels::Vector{String}
)::TTestResult
    m1 = mean(x1)
    m2 = mean(x2)
    ms = [m1, m2]
    v1 = var(x1)
    v2 = var(x2)
    sds = sqrt.([v1, v2])
    n1 = length(x1)
    n2 = length(x2)
    ns = [n1, n2]

    Δ = m1 - m2
    ν = n1 + n2 - 2
    σ = sqrt(((n1-1)*v1 + (n2-1)*v2) / ν)
    t_stat = Δ / (σ * sqrt(1/n1 + 1/n2))
    pval = 2.0 * D.ccdf(TDist(ν), abs(t_stat))
    distr_t = TDist.(ns .- 1)
    qm = [quantile(t, p) for (t, p) in Iterators.product(distr_t, [0.025, 0.975])]
    ci = ms .+ qm .* sds ./ sqrt.(ns)
    cohen_d = abs(Δ) / sqrt((v1 + v2) / 2.0)

    res = TTestResult(
        :simple, levels, ns, ms, sds, ci,
        t_stat, ν, pval, cohen_d
    )
    return res
end

function t_test_welch(
    x1::FloatVector,
    x2::FloatVector,
    levels::Vector{String}
)::TTestResult
    m1 = mean(x1)
    m2 = mean(x2)
    ms = [m1, m2]
    v1 = var(x1)
    v2 = var(x2)
    sds = sqrt.([v1, v2])
    n1 = length(x1)
    n2 = length(x2)
    ns = [n1, n2]
    vn1 = v1 / n1
    vn2 = v2 / n2
    
    Δ = m1 - m2
    σ = sqrt(vn1 + vn2)
    t_stat = Δ / σ
    ν = (vn1 + vn2)^2 / ((vn1^2/(n1-1)) + (vn2^2/(n2-1)))
    pval = 2.0 * D.ccdf(TDist(ν), abs(t_stat))
    distr_t = TDist.(ns .- 1)
    qm = [quantile(t, p) for (t, p) in Iterators.product(distr_t, [0.025, 0.975])]
    ci = ms .+ qm .* sds ./ sqrt.(ns)
    cohen_d = abs(Δ) / sqrt((v1 + v2) / 2.0)
    
    res = TTestResult(:welch, levels, ns, ms, sds, ci,
                      t_stat, ν, pval, cohen_d)
    return res
end

"""
    TTestModel

A `TestModel` type for t-test

# Members

- `X`: model matrix
- `y`: vector to compare its means
- `type`: either :one, :paired, :simple, or :welch
- `levels`: `String` vector of the unique levels (only 1 or 2 elements are allowed)
"""
struct TTestModel{T<:Real} <: TestModel
    x1::AbstractVector{T}
    x2::AbstractVector{T}
    type::Symbol
    levels::Vector{String}
    function TTestModel(
        X::Matrix{T},
        y::AbstractVector{T},
        type::Symbol,
        levels::Vector{String};
        μ0::T = 0.0
    ) where T
        if type == :one
            new{T}(y, [μ0], :one, levels)
        else
            dummy = X[:,2]
            x1 = y[dummy .== 0]
            x2 = y[dummy .== 1]
            new{T}(x1, x2, type, levels)
        end
    end
end

function fit(obj::TTestModel)
    if obj.type == :one
        t_test(obj.x1, μ0=obj.x2[1])
    elseif obj.type == :paired
        t_test(obj.x1, obj.x2, paired=true, levels=obj.levels)
    else
        t_test(obj.x1, obj.x2, paired=false, varequal=(obj.type == :simple), levels=obj.levels)
    end
end

function fit(
        ::Type{TTestModel},
        X::AbstractMatrix{<:Real},
        y::AbstractVector{<:Real},
        type::Symbol,
        levels::Vector{String};
        μ0::Real = 0.0
    )
    ttest_obj = TTestModel(X, y, type, levels, μ0=μ0)
    fit(ttest_obj)
end
function fit(
        ::Type{TTestModel},
        X::AbstractMatrix{<:Real},
        y::AbstractVector{<:Real},
        type::Symbol
    )
    @info "Using the first and second sample names as 'Sample1' and 'Sample2'"
    fit(TTestModel, X, y, type, ["Sample1", "Sample2"], μ0=0.0)
end

"""
    t_test(X, y, paired::Bool=false, varequal::Bool=false)

An alias for `fit(TTestModel, X, y, levels)`

The arguments `X` and `y` can be a `Matrix` and a `Vector` or a `Formula` and `DataFrame`.

The argument `type` must be either `:one`, `:paired`, `:simple`, or `welch`.

The argument `levels` can be a `Vector{String}`.

The keyword argument `μ0` can be a `Real` to use one sample t-test (default = 0.0). Ignored if the `type` is not `:one`.
"""
function t_test(
    X::AbstractMatrix{<:Real},
    y::AbstractVector{<:Real},
    paired::Bool=false,
    varequal::Bool=false
)
    if paired
        fit(TTestModel, X, y, :paired)
    elseif varequal
        fit(TTestModel, X, y, :simple)
    else
        fit(TTestModel, X, y, :welch)
    end
end
function t_test(
    formula::FormulaTerm{Term,Term},
    df::DataFrame,
    paired::Bool=false,
    varequal::Bool=false
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

    # perform t-test
    if paired
        fit(TTestModel, X, y, :paired, levels)
    elseif varequal
        fit(TTestModel, X, y, :simple, levels)
    else
        fit(TTestModel, X, y, :welch, levels)
    end
end
