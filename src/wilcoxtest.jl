mutable struct WilcoxTestResult<:TestResult
    type::Symbol  # either :paired or :independent
    is_corrected::Bool
    levels::Union{Vector{String}, Nothing}
    n_obs::IntVector
    means::FloatVector
    sds::FloatVector
    statistic::Float64
    pval::Float64
    cohen_d::Union{Float64, Nothing}
end

function show(io::IO, wtr::WilcoxTestResult)
    if wtr.type == :paired
        if wtr.is_corrected
            println(
                io,
                "\n        Wilcoxon signed rank sum test with continuity correction\n"
            )
        else
            println(io, "\n        Wilcoxon signed rank sum exact test\n")
        end
        @printf(io, "Difference (mean ± sd): %.2f ± %.2f\n",
                wtr.means[1], wtr.sds[1])
    elseif wtr.type == :independent
        if wtr.is_corrected
            println(
                io,
                "\n        Wilcoxon rank sum test with continuity correction\n"
            )
        else
            println(io, "\n        Wilcoxon rank sum exact test\n")
        end
        println(io, "Data (mean ± sd)")
        @printf(io, "  %s: %.2f ± %.2f\n",
                wtr.levels[1], wtr.means[1], wtr.sds[1])
        @printf(io, "  %s: %.2f ± %.2f\n",
                wtr.levels[2], wtr.means[2], wtr.sds[2])
    end
    println(io, "\nResult:\n", coeftable(wtr))
end

function coeftable(wtr::WilcoxTestResult)
    if (wtr.type == :paired)
        CoefTable(
            [wtr.statistic, wtr.pval],
            ["W", "P-value"],
            [""],
            2, 1
        )
    elseif (wtr.type == :independent)
        CoefTable(
            [wtr.statistic, wtr.pval, wtr.cohen_d],
            ["W", "P-value", "Cohen's d"],
            [""],
            2, 1
        )
    else
        error("invalid Wilcoxon test type $wtr.type")
    end
end


"""
        WilcoxDist(n1, n2)

Distribution of the Wilcox statistic
"""
struct WilcoxDist <: Distributions.DiscreteUnivariateDistribution
    n1::Int
    n2::Int
    w::Array{Float64, 3}
end
Distributions.@distr_support WilcoxDist 0 n1*n2
WilcoxDist(n1::Int, n2::Int) = WilcoxDist(
    n1,
    n2,
    fill(-1.0, (max(n1, n2, 50), 
                max(n1, n2, 50),
                ceil(Int, max(n1,50) * max(n2,50) / 2)))
)
Base.show(io::IO, wd::WilcoxDist) = print(io, "Wilcox distribution with (", wd.n1, ", ", wd.n2, ")\n")

# NOTE: Julia starts from 1, not from 0
function cwilcox(w::Array{Float64, 3}, k::Int, m::Int, n::Int)::Float64
    mn = m * n
    if (k < 0) | (k > mn+1)
        return 0.0
    end
    c = floor(Int, mn / 2)
    if k > c
        k = mn - k  # hence k <= floor(mn/2)
    end
    i = min(m, n)
    j = max(m, n)  # hence i <= j
    if j == 0  # hence i == 0
        if k == 0
            return 1.0
        else
            return 0.0
        end
    end

    if (j > 0) & (k < j)
        return cwilcox(w, k, i, k)
    end
    if w[i,j,k] < 0
        w[i,j,k] = cwilcox(w, k-j, i-1, j) + cwilcox(w, k, i, j-1)
    end
    # @printf(stdout, "(i,j,k) = (%d,%d,%d) | c = %.2f\n", i, j, k, w[i,j,k])
    return w[i,j,k]
end

function cdf(d::WilcoxDist, x::Float64)
    n12 = d.n1 * d.n2
    c = binomial(big(d.n1 + d.n2), big(d.n2))
    if x < 0.0
        return 0.0
    elseif x ≥ n12
        return 1.0
    else
        p = 0.0
        if x ≤ (n12/2.0)
            for i = 0:floor(Int, x)
                p += cwilcox(d.w, i, d.n1, d.n2) / c
            end
        else
            x = n12 - x
            for i = 0:(floor(Int, x)-1)
                p += cwilcox(d.w, i, d.n1, d.n2) / c
            end
        end
        return p
    end
end

ccdf(d::WilcoxDist, x::Float64) = 1.0 - cdf(d, x)


"""
        WilcoxSignedDist(n1, n2)

Distribution of the Wilcox statistic
"""
struct WilcoxSignedDist <: Distributions.DiscreteUnivariateDistribution
    n::Int
    w::Array{Float64}
end
Distributions.@distr_support WilcoxSignedDist 0 n
WilcoxSignedDist(n::Int) = WilcoxSignedDist(
    n,
    zeros(ceil(Int, max(n,50)*(max(n,50)+1)/2))
)
Base.show(io::IO, wsd::WilcoxSignedDist) = print(io, "Wilcox signed distribution with ", wsd.n, "\n")

function cwilcoxsigned(w::Array{Float64}, k::Int, n::Int)::Float64
    u = n * (n+1) / 2.0
    c = floor(Int, u/2)

    if (k < 1) | (k > u)
        return 0.0
    end

    if k > c
        k = u - k
    end

    if n == 1
        return 1.0
    end

    if w[1] == 1.0
        return w[k]
    else
        w[1] = 1.0
        w[2] = 1.0
        for j = 3:(n+2)
            i = min(floor(Int, j*(j+1)/2), c)
            while i ≥ j
                w[i] += w[i-j+1]
                i -= 1
            end
        end
        return w[k]
    end
end

function cdf(d::WilcoxSignedDist, x::Float64)
    if x < 0.0
        return 0.0
    elseif x > (d.n*(d.n+1)/2.0)
        return 1.0
    else
        nf = float(d.n)
        f = exp(- nf * log(2.0))
        p = 0.0
        if x ≤ (d.n*(d.n+1)/4.0)
            for i = 1:(ceil(Int, x)+1)
                p += cwilcoxsigned(d.w, i, d.n) * f
            end
        else
            x = d.n*(d.n+1)/2.0 - x
            for i = 1:ceil(Int, x)
                p += cwilcoxsigned(d.w, i, d.n) * f
            end
        end
        return p
    end
end

ccdf(d::WilcoxSignedDist, x::Float64) = 1.0 - cdf(d, x)



"""
    wilcox_test(x1, x2, paired::Bool=false)

Wilcoxon (signed) rank sum test

The arguments `x1` and `x2` can be numeric `Vector`.

Wilcoxon signed rank sum test will be performed if the keyword argument `paired` is `true`. Otherwise, Wilcoxon rank sum test will be performed.
"""
function wilcox_test(
    x1::NumVector,
    x2::NumVector;
    levels::Union{Vector{String}, Nothing}=nothing,
    paired::Bool=false,
    exact::Union{Bool, Nothing}=nothing,
    μ0::AbstractFloat=0.0
)::WilcoxTestResult
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

    # paired or independent
    if paired
        n1 == n2 || error("size of $levels[1] and $levels[2] must be the same: ($n1, $n2)")
        if isnothing(exact)
            exact = n1 < 50
        end
        wilcox_test_signed(x1, x2, exact, μ0)
    else
        if isnothing(exact)
            exact = (n1 < 50) & (n2 < 50)
        end
        wilcox_test_independent(x1, x2, levels, exact, μ0)
    end
end

function wilcox_test_signed(
    x1::FloatVector,
    x2::FloatVector,
    exact::Bool,
    μ0::AbstractFloat
)::WilcoxTestResult
    n = length(x1)
    x = x1 .- x2 .- μ0

    x = filter(a -> a ≠ 0.0, x)
    zeroes = length(x) < n
    n = length(x)  # redefinition

    r = tiedrank(abs.(x))
    w_stat = sum(r[x .> 0])

    n_ties = values(countmap(r))
    ties = false
    if length(n_ties) < n
        ties = true
    end

    pval = 1.0
    is_corrected = false
    if exact & !ties & !zeroes
        p = 1.0
        d = WilcoxSignedDist(n)
        if w_stat > (n*(n+1)/4.0)
            p = cdf(d, w_stat-1.0)
        else
            p = ccdf(d, w_stat)
        end
        pval = min(2.0 * p, 1.0)
    else
        is_corrected = true
        z = w_stat - n*(n+1)/4.0
        σ = sqrt(n*(n+1)*(2*n+1)/24.0 -
            sum(n_ties.^3 .- n_ties)/48.0)
        correction = 0.5 * sign(z)
        z = (z - correction) / σ
        norm_std = Normal()
        pval = 2.0 * min(
            Distributions.cdf(norm_std, z),
            Distributions.ccdf(norm_std, z)
        )
    end

    res = WilcoxTestResult(
        :paired, is_corrected, nothing, [n], [mean(x)],
        [std(x)], w_stat, pval, nothing
    )
    return res
end

function wilcox_test_independent(
    x1::FloatVector,
    x2::FloatVector,
    levels::Vector{String},
    exact::Bool,
    μ0::AbstractFloat
)::WilcoxTestResult
    n1 = length(x1)
    n2 = length(x2)
    r = tiedrank(vcat(x1 .- μ0, x2))
    n_ties = values(countmap(r))
    
    w_stat = sum(r[1:n1]) - n1 * (n1 + 1) / 2
    ties = false
    if length(n_ties) < (n1+n2)
        ties = true
    end

    pval = 1.0
    is_corrected = false
    if exact & !ties
        p = 1.0
        d = WilcoxDist(n1, n2)
        if w_stat > (n1 * n2 / 2)
            p = cdf(d, w_stat-1.0)
        else
            p = ccdf(d, w_stat)
        end
        pval = min(2.0*p, 1.0)
    else
        is_corrected = true
        z = w_stat - n1 * n2 / 2.0
        σ = sqrt((n1*n2/12.0) * 
            ((n1+n2+1.0) - sum(n_ties.^3 .- n_ties) / ((n1+n2) * (n1+n2-1.0))))
        correction = 0.5 * sign(z)
        z = (z - correction) / σ
        norm_std = Normal()
        pval = 2.0 * min(
            Distributions.cdf(norm_std, z),
            Distributions.ccdf(norm_std, z)
        )
    end

    m1 = mean(x1)
    m2 = mean(x2)
    v1 = var(x1)
    v2 = var(x2)
    cohen_d = abs(m1-m2) / sqrt((v1 + v2) / 2.0)
    
    res = WilcoxTestResult(
        :independent, is_corrected, levels, [n1, n2], [m1, m2],
        [sqrt(v1), sqrt(v2)], w_stat, pval, cohen_d
    )
    return res
end

"""
    WilcoxTestModel

A `WilcoxestModel` type for Wilcoxon (signed) rank sum test

# Members

- `X`: model matrix
- `y`: vector to compare its means
- `type`: either `:paired` or `:independent`
- `levels`: `String` vector of the unique levels (only 1 or 2 elements are allowed)
"""
struct WilcoxTestModel{T<:Real} <: TestModel
    x1::AbstractVector{T}
    x2::AbstractVector{T}
    type::Symbol
    levels::Vector{String}
    function WilcoxTestModel(
        X::Matrix{T},
        y::AbstractVector{T},
        type::Symbol,
        levels::Vector{String};
        μ0::T = 0.0
    ) where T
        dummy = X[:,2]
        x1 = y[dummy .== 0]
        x2 = y[dummy .== 1]
        new{T}(x1, x2, type, levels)
    end
end

function fit(obj::WilcoxTestModel)
    if obj.type == :paired
        wilcox_test(obj.x1, obj.x2, paired=true, levels=obj.levels)
    else
        wilcox_test(obj.x1, obj.x2, paired=false, levels=obj.levels)
    end
end

function fit(
        ::Type{WilcoxTestModel},
        X::AbstractMatrix{<:Real},
        y::AbstractVector{<:Real},
        type::Symbol,
        levels::Vector{String};
        μ0::Real = 0.0
    )
    wilcox_test_obj = WilcoxTestModel(X, y, type, levels, μ0=μ0)
    fit(wilcox_test_obj)
end
function fit(
        ::Type{WilcoxTestModel},
        X::AbstractMatrix{<:Real},
        y::AbstractVector{<:Real},
        type::Symbol
    )
    @info "Using the first and second sample names as 'Sample1' and 'Sample2'"
    fit(WilcoxTestModel, X, y, type, ["Sample1", "Sample2"], μ0=0.0)
end

"""
    wilcox_test(X, y, paired::Bool=false)

An alias for `fit(WilcoxTestModel, X, y, levels)`

The arguments `X` and `y` can be a `Matrix` and a `Vector` or a `Formula` and `DataFrame`.

The argument `type` must be either `:paired` or `:independent`.

The argument `levels` can be a `Vector{String}`.
"""
function wilcox_test(
    X::AbstractMatrix{<:Real},
    y::AbstractVector{<:Real},
    paired::Bool=false
)
    if paired
        fit(WilcoxTestModel, X, y, :paired)
    else
        fit(WilcoxTestModel, X, y, :independent)
    end
end
function wilcox_test(
    formula::FormulaTerm{Term,Term},
    df::DataFrame,
    paired::Bool=false
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
        fit(WilcoxTestModel, X, y, :paired, levels)
    else
        fit(WilcoxTestModel, X, y, :independent, levels)
    end
end
