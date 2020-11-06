"""
    t_test(x, mu::AbstractFloat=0.0)

One sample t-test

The argument `x` can be numeric `Vector`.

One sample t-test will be performed to compare between the mean of `x1` and `mu`.
"""
function t_test(x::NumVector; mu::AbstractFloat=0.0)::Dict
    # check types
    if !(x isa FloatVector)
        @warn "x1 is converted to Float64"
        x = convert(Vector{Float64}, x)
    end
    
    # check length
    n = length(x)
    n > 1 || error("size of x is not enough: actual size is ", n)

    m = mean(x)
    s = std(x)
    t_stat = abs(m - mu) / (s / sqrt(n))
    nu = n - 1
    pval = 2.0 * ccdf(TDist(nu), t_stat)

    res = Dict(
        "statistic" => t_stat,
        "df" => nu,
        "P-value" => pval,
    )
    println("\n        One sample t-test\n")
    println("Data (mean ± sd)")
    @printf("  x: %.2f ± %.2f\n", m, s)
    @printf("  difference from %.2f\n", mu)
    println("\nResult")
    @printf(
        "  t = %.3f, df = %.2f, P-value = %.3f\n",
        t_stat, nu, pval
    )
    print("\n")
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
)::Dict
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
    n1 > 1 || error("size of x1 is not enough: actual size is ", n1)
    n2 > 1 || error("size of x2 is not enough: actual size is ", n2)

    # paired, simple, or Welch
    if paired
        n1 == n2 || error("size of x1 and x2 must be the same")
        t_test_paired(x1, x2)
    elseif varequal
        t_test_var_equal(x1, x2)
    else
        t_test_welch(x1, x2)
    end
end

function t_test_paired(x1::FloatVector, x2::FloatVector)::Dict
    n = length(x1)
    diff = x1 .- x2
    mean_diff = abs(mean(diff))
    sd_diff = std(diff)
    t_stat = mean_diff / (sd_diff / sqrt(n))
    nu = n - 1
    pval = 2.0 * ccdf(TDist(nu), t_stat)

    res = Dict(
        "statistic" => t_stat,
        "df" => nu,
        "P-value" => pval,
    )
    println("\n        Paired t-test\n")
    @printf("Difference (mean ± sd): %.2f ± %.2f\n", mean_diff, sd_diff)
    println("\nResult")
    @printf(
        "  t = %.3f, df = %.2f, P-value = %.3f\n",
        t_stat, nu, pval
    )
    print("\n")
    return res
end

function t_test_var_equal(x1::FloatVector, x2::FloatVector)::Dict
    m1 = mean(x1)
    m2 = mean(x2)
    v1 = var(x1)
    v2 = var(x2)
    n1 = length(x1)
    n2 = length(x2)

    diff_means = abs(m1 - m2)
    nu = n1 + n2 - 2
    sd_hat = ((n1-1)*v1 + (n2-1)*v2) / nu
    t_stat = diff_means / (sd_hat * sqrt(1/n1 + 1/n2))
    pval = 2.0 * ccdf(TDist(nu), t_stat)
    cohen_d = diff_means / sqrt((v1 + v2) / 2.0)

    res = Dict(
        "statistic" => t_stat,
        "df" => nu,
        "P-value" => pval,
        "Cohen's d" => cohen_d,
    )
    println("\n        Simple Two Sample t-test\n")
    println("Data (mean ± sd)")
    @printf("  x1: %.2f ± %.2f\n", m1, sqrt(v1))
    @printf("  x2: %.2f ± %.2f\n", m2, sqrt(v2))
    println("\nResult")
    @printf(
        "  t = %.3f, df = %.2f, P-value = %.3f, Cohen's d = %.2f\n",
        t_stat, nu, pval, cohen_d
    )
    print("\n")
    return res
end

function t_test_welch(x1::FloatVector, x2::FloatVector)::Dict
    m1 = mean(x1)
    m2 = mean(x2)
    v1 = var(x1)
    v2 = var(x2)
    n1 = length(x1)
    n2 = length(x2)
    vn1 = v1 / n1
    vn2 = v2 / n2
    
    diff_means = abs(m1 - m2)
    sd_hat = sqrt(vn1 + vn2)
    t_stat = diff_means / sd_hat
    nu = (vn1 + vn2)^2 / ((vn1^2/(n1-1)) + (vn2^2/(n2-1)))
    pval = 2.0 * ccdf(TDist(nu), t_stat)
    cohen_d = diff_means / sqrt((v1 + v2) / 2.0)
    
    res = Dict(
        "statistic" => t_stat,
        "df" => nu,
        "P-value" => pval,
        "Cohen's d" => cohen_d,
    )
    println("\n        Welch Two Sample t-test\n")
    println("Data (mean ± sd)")
    @printf("  x1: %.2f ± %.2f\n", m1, sqrt(v1))
    @printf("  x2: %.2f ± %.2f\n", m2, sqrt(v2))
    println("\nResult")
    @printf(
        "  t = %.3f, df = %.2f, P-value = %.3f, Cohen's d = %.2f\n",
        t_stat, nu, pval, cohen_d
    )
    print("\n")
    return res
end
