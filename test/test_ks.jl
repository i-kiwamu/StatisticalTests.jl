@testset "Kolmogorov-Smirnov test" begin
    x = [
        392,585,410,414,384,423,456,549,455,406,470,402,416,431,435,413,468,435,
        432,438,425,392,410,498,475,431,460,440,420,422,409,417,370,500,501,478,
        454,444,447,393,446,454,500,433,444,427,412,433,460,434
    ]
    res_kstest = ks_test(x, Normal(mean(x), std(x)))
    @test res_kstest.statistic ≈ 0.119089045870544
    @test res_kstest.n_obs == 50
    @test res_kstest.pval ≈ 0.477419559223632
end
