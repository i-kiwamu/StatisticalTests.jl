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

function show(io::IO, ttr::TTestResult)
    if ttr.type == :one
        println("\n        One sample t-test\n")
        println("Data (mean ± sd)")
        @printf("  x: %.2f ± %.2f\n", ttr.means[1], ttr.sds[1])
        @printf("  difference from %.2f\n", ttr.means[2])
        println("\nResult")
        @printf(
            "  t = %.3f, df = %d, P-value = %.3f\n",
            ttr.statistic, ttr.df, ttr.pval
        )
        print("\n")
    elseif ttr.type == :paired
        println("\n        Paired t-test\n")
        @printf("Difference (mean ± sd): %.2f ± %.2f\n",
                ttr.means[1], ttr.sds[1])
        println("\nResult")
        @printf(
            "  t = %.3f, df = %d, P-value = %.3f\n",
            ttr.statistic, ttr.df, ttr.pval
        )
        print("\n")
    elseif ttr.type == :simple
        println("\n        Simple Two Sample t-test\n")
        println("Data (mean ± sd)")
        @printf("  x1: %.2f ± %.2f\n", ttr.means[1], ttr.sds[1])
        @printf("  x2: %.2f ± %.2f\n", ttr.means[2], ttr.sds[2])
        println("\nResult")
        @printf(
            "  t = %.3f, df = %d, P-value = %.3f, Cohen's d = %.2f\n",
            ttr.statistic, ttr.df, ttr.pval, ttr.cohen_d
        )
        print("\n")
    elseif ttr.type == :welch
        println("\n        Welch Two Sample t-test\n")
        println("Data (mean ± sd)")
        @printf("  x1: %.2f ± %.2f\n", ttr.means[1], ttr.sds[1])
        @printf("  x2: %.2f ± %.2f\n", ttr.means[2], ttr.sds[2])
        println("\nResult")
        @printf(
            "  t = %.3f, df = %.2f, P-value = %.3f, Cohen's d = %.2f\n",
            ttr.statistic, ttr.df, ttr.pval, ttr.cohen_d
        )
        print("\n")    
    else
        error("invalid t-test type $ttr.type")
    end
end

"""
    t_test(x, mu::AbstractFloat=0.0)

One sample t-test

The argument `x` can be numeric `Vector`.

One sample t-test will be performed to compare between the mean of `x1` and `mu`.
"""
function t_test(x::NumVector; mu::AbstractFloat=0.0)::TTestResult
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
    t_stat = (μ - mu) / (σ / sqrt(n))
    ν = n - 1
    pval = 2.0 * ccdf(TDist(ν), abs(t_stat))

    res = TTestResult(:one, [n], [μ, mu], [σ], t_stat, ν, pval, nothing)
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
