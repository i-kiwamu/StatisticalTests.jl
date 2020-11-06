using Test, RDatasets
using StatisticalTests

@testset "Welch two sample t-test" begin
    df_sleep = dataset("datasets", "sleep")
    res_ttest = t_test(
        df_sleep.Extra[df_sleep.Group .== "1"],
        df_sleep.Extra[df_sleep.Group .== "2"]
    )
    test_show(res_ttest)
    @test isapprox(res_ttest["statistic"], 1.8608)
    @test isapprox(res_ttest["df"], 17.776)
    @test isapprox(res_ttest["P-value"], 0.079394)
end
