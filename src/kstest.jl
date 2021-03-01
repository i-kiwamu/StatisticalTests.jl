mutable struct KSTestResult <: TestResult
    n_obs::Int
    statistic::Float64
    pval::Float64
end

function show(io::IO, kstr::KSTestResult)
    println(io, "\n        One Sample Kolmogorov-Smirnov test")
    println(io, "\n Result:\n", coeftable(kstr))
end

function coeftable(kstr::KSTestResult)
    CoefTable(
        [kstr.n_obs, kstr.statistic, kstr.pval],
        ["N", "D", "P-value"],
        [""],
        3, 2
    )
end

function pkstwo(x::Float64, tol = 1e-06)
    k_max = floor(Int, sqrt(2.0 - log(tol)))
    if x < 1
        z = - (π^2 / 8) / (x^2)
        w = log(x)
        s = 0.0
        k = 1
        while k < k_max
            s += exp(k^2 * z - w)
            k += 2
        end
        return s * sqrt(2 * π)
    else
        z = -2 * x^2
        s = -1.0
        old = 0.0
        new = 1.0
        while abs(old - new) > tol
            old = new
            new += 2 * s * exp(z * k^2)
            s *= -1.0
            k += 1
        end
        return new
    end
end

function ccdfNonExact(d::KSDist, x::Float64)
    return 1 - pkstwo(sqrt(d.n) * x)
end

"""
    ks_test(x, distr)

One sample Kolmogorov-Smirnov test

The argument `x` can be a numeric `Vector`.

The argument `distr` can be a `Distributions.UnivariateDistribution`.
"""
function ks_test(
    x::NumVector,
    distr::UnivariateDistribution;
    exact::Union{Bool, Nothing}=nothing
)::KSTestResult
    n = length(x)
    ties = false
    if length(unique(x)) < n
        @warn "ties should not be present for the Kolmogorov-Smirnov test"
        ties = true
    end

    if isnothing(exact)
        exact = (n < 100) & !ties
    end

    x_density_sorted = cdf.(distr, sort(x)) - collect(0:(n-1))/n
    statistic = max(maximum(x_density_sorted), maximum(1/n .- x_density_sorted))

    if exact
        pval = ccdf(KSDist(n), statistic)
    else
        pval = ccdfNonExact(KSDist(n), statistic)
    end

    res = KSTestResult(n, statistic, pval)
    return res
end
