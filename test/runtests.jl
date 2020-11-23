using StatisticalTests
using Test, RDatasets, StatsModels

test_show(x) = show(IOBuffer(), x)

@testset "One sample t-test" begin
    x = [2.55, 1.94, 2.52, 2.47, 2.00, 3.04, 2.82, 3.45, 3.82, 0.89]
    res_ttest = t_test(x, μ0=0.0)
    test_show(res_ttest)
    @test isapprox(res_ttest.statistic, 9.73641835357921)
    @test res_ttest.dof == 9
    @test isapprox(res_ttest.pval, 4.46718590474277e-06)
end

@testset "Paired t-test" begin
    x1 = [-2.0, -3.2, 1.2, 5.8, 19.6, 19.9, 28.9, 31.1, 27.9, 21.9, 8.8, 7.5]
    x2 = [-1.8, -2.9, 0.9, 1.5, 11.5, 18.8, 23.9, 24.9, 23.8, 19.4, 8.8, 4.8]
    res_ttest = t_test(x1, x2, paired=true)
    # test_show(res_ttest)
    @test isapprox(res_ttest.statistic, 3.52842907954346)
    @test res_ttest.dof == 11
    @test isapprox(res_ttest.pval, 0.0047278741207561)
end

@testset "Simple two sample t-test" begin
    x1 = [21.8, 21.7, 32.0, 30.2, 27.4]
    x2 = [12.3, 19.7, 21.9, 23.6]
    res_ttest = t_test(x1, x2, varequal=true)
    # test_show(res_ttest)
    @test isapprox(res_ttest.statistic, 2.22997773895153)
    @test res_ttest.dof == 7
    @test isapprox(res_ttest.pval, 0.0609728704866215)
end

@testset "Welch two sample t-test" begin
    df_sleep = dataset("datasets", "sleep")
    res_ttest = t_test(@formula(Extra ~ Group), df_sleep)
    test_show(res_ttest)
    @test isapprox(res_ttest.statistic, -1.86081346748685)
    @test isapprox(res_ttest.dof, 17.7764735161785)
    @test isapprox(res_ttest.pval, 0.0793941401873582)
end

@testset "F-test 1" begin
    x1 = [19.4, 23.9, 21.0, 22.8, 15.8]
    x2 = [34.9, 27.1, 6.9, 13.2, 7.3]
    res_ftest = f_test(x1, x2)
    test_show(res_ftest)
    @test isapprox(res_ftest.statistic, 15.5838287752675)
    @test res_ftest.dofs[1] == 4
    @test res_ftest.dofs[2] == 4
    @test isapprox(res_ftest.pval, 0.0209393152025605)
end

@testset "F-test 2" begin
    df_sleep = dataset("datasets", "sleep")
    res_ftest = f_test(@formula(Extra ~ Group), df_sleep)
    test_show(res_ftest)
    @test isapprox(res_ftest.statistic, 1.2525950355841)
    @test res_ftest.dofs[1] == 9
    @test res_ftest.dofs[2] == 9
    @test isapprox(res_ftest.pval, 0.742719931726045)
end
