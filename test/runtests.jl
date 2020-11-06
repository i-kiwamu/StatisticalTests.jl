using Test, RDatasets
using StatisticalTests

test_show(x) = show(IOBuffer(), x)

@testset "Welch two sample t-test" begin
    df_sleep = dataset("datasets", "sleep")
    res_ttest = t_test(
        df_sleep.Extra[df_sleep.Group .== "1"],
        df_sleep.Extra[df_sleep.Group .== "2"]
    )
    test_show(res_ttest)
    @test isapprox(res_ttest["statistic"], 1.86081346748685)
    @test isapprox(res_ttest["df"], 17.7764735161785)
    @test isapprox(res_ttest["P-value"], 0.0793941401873582)
end
