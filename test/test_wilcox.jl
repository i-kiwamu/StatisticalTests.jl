@testset "Wilcoxon rank sum exact test" begin
    x1 = [0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46]
    x2 = [1.15, 0.88, 0.90, 0.74, 1.21]
    res_wtest = wilcox_test(x1, x2)
    test_show(res_wtest)
    @test res_wtest.statistic ≈ 35
    @test res_wtest.pval ≈ 0.254412254412254
end

@testset "Wilcoxon rank sum test with continuity correction" begin
    df_sleep = dataset("datasets", "sleep")
    res_wtest = wilcox_test(@formula(Extra ~ Group), df_sleep)
    test_show(res_wtest)
    @test res_wtest.statistic ≈ 25.5
    @test res_wtest.pval ≈ 0.0693275754336266
end

@testset "Wilcoxon signed rank sum exact test" begin
    x1 = [1.83,  0.50,  1.62,  2.48, 1.68, 1.88, 1.55, 3.06, 1.30]
    x2 = [0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29]
    res_wtest = wilcox_test(x1, x2, paired=true)
    test_show(res_wtest)
    @test res_wtest.statistic ≈ 40
    @test res_wtest.pval ≈ 0.0390625
end

@testset "Wilcoxon signed rank sum test with continuity correction" begin
    x1 = [-2.0, -3.2, 1.2, 5.8, 19.6, 19.9, 28.9, 31.1, 27.9, 21.9, 8.8, 7.5]
    x2 = [-1.8, -2.9, 0.9, 1.5, 11.5, 18.8, 23.9, 24.9, 23.8, 19.4, 8.8, 4.8]
    res_wtest = wilcox_test(x1, x2, paired=true)
    test_show(res_wtest)
    @test res_wtest.statistic ≈ 62
    @test res_wtest.pval ≈ 0.011278190116074
end
