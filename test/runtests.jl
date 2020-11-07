using StatisticalTests
using Test, RDatasets

test_show(x) = show(IOBuffer(), x)

@testset "Paired t-test" begin
    x1 = [-2.0, -3.2, 1.2, 5.8, 19.6, 19.9, 28.9, 31.1, 27.9, 21.9, 8.8, 7.5]
    x2 = [-1.8, -2.9, 0.9, 1.5, 11.5, 18.8, 23.9, 24.9, 23.8, 19.4, 8.8, 4.8]
    res_ttest = t_test(x1, x2, paired=true)
    test_show(res_ttest)
    @test isapprox(res_ttest.statistic, 3.52842907954346)
    @test res_ttest.df == 11
    @test isapprox(res_ttest.pval, 0.0047278741207561)
end

@testset "Simple two sample t-test" begin
    x1 = [21.8, 21.7, 32.0, 30.2, 27.4]
    x2 = [12.3, 19.7, 21.9, 23.6]
    res_ttest = t_test(x1, x2, varequal=true)
    test_show(res_ttest)
    @test isapprox(res_ttest.statistic, 2.22997773895153)
    @test res_ttest.df == 7
    @test isapprox(res_ttest.pval, 0.0609728704866215)
end

@testset "Welch two sample t-test" begin
    df_sleep = dataset("datasets", "sleep")
    res_ttest = t_test(
        df_sleep.Extra[df_sleep.Group .== "1"],
        df_sleep.Extra[df_sleep.Group .== "2"]
    )
    test_show(res_ttest)
    @test isapprox(res_ttest.statistic, -1.86081346748685)
    @test isapprox(res_ttest.df, 17.7764735161785)
    @test isapprox(res_ttest.pval, 0.0793941401873582)
end
