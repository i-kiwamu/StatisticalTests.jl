@testset "F-test 1" begin
    x1 = [19.4, 23.9, 21.0, 22.8, 15.8]
    x2 = [34.9, 27.1, 6.9, 13.2, 7.3]
    res_ftest = f_test(x1, x2)
    test_show(res_ftest)
    @test res_ftest.statistic ≈ 15.5838287752675
    @test res_ftest.dofs[1] == 4
    @test res_ftest.dofs[2] == 4
    @test res_ftest.pval ≈ 0.0209393152025605
end

@testset "F-test 2" begin
    df_sleep = dataset("datasets", "sleep")
    res_ftest = f_test(@formula(Extra ~ Group), df_sleep)
    test_show(res_ftest)
    @test res_ftest.statistic ≈ 1.2525950355841
    @test res_ftest.dofs[1] == 9
    @test res_ftest.dofs[2] == 9
    @test res_ftest.pval ≈ 0.742719931726045
end
