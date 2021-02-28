var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = StatisticalTests","category":"page"},{"location":"#StatisticalTests","page":"Home","title":"StatisticalTests","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [StatisticalTests]","category":"page"},{"location":"#StatisticalTests.FTestModel","page":"Home","title":"StatisticalTests.FTestModel","text":"FTestModel\n\nA TestModel type for f-test\n\nMembers\n\nX: model matrix\ny: vector to compare its means\nlevels: String vector of the unique levels (only 2 elements are allowed)\n\n\n\n\n\n","category":"type"},{"location":"#StatisticalTests.TTestModel","page":"Home","title":"StatisticalTests.TTestModel","text":"TTestModel\n\nA TestModel type for t-test\n\nMembers\n\nX: model matrix\ny: vector to compare its means\ntype: either :one, :paired, :simple, or :welch\nlevels: String vector of the unique levels (only 1 or 2 elements are allowed)\n\n\n\n\n\n","category":"type"},{"location":"#StatisticalTests.f_test-Tuple{AbstractArray{T,1} where T<:Number,AbstractArray{T,1} where T<:Number}","page":"Home","title":"StatisticalTests.f_test","text":"f_test(x1, x2)\n\nF-test\n\nThe arguments x1 and x2 can be numeric Vector.\n\n\n\n\n\n","category":"method"},{"location":"#StatisticalTests.f_test-Tuple{AbstractArray{var\"#s141\",2} where var\"#s141\"<:Real,AbstractArray{var\"#s140\",1} where var\"#s140\"<:Real}","page":"Home","title":"StatisticalTests.f_test","text":"f_test(X, y)\n\nAn alias for fit(FTestModel, X, y, levels)\n\nThe arguments X and y can be a Matrix and a Vector or a Formula and DataFrame.\n\nThe argument levels can be a Vector{String}.\n\n\n\n\n\n","category":"method"},{"location":"#StatisticalTests.isna-Tuple{Any}","page":"Home","title":"StatisticalTests.isna","text":"isna(x)\n\nIndicate whether x is missing or NaN.\n\n\n\n\n\n","category":"method"},{"location":"#StatisticalTests.ks_test-Tuple{AbstractArray{T,1} where T<:Number,Distributions.Distribution{Distributions.Univariate,S} where S<:Distributions.ValueSupport}","page":"Home","title":"StatisticalTests.ks_test","text":"ks_test(x, distr)\n\nOne sample Kolmogorov-Smirnov test\n\nThe argument x can be a numeric Vector.\n\nThe argument distr can be a Distributions.UnivariateDistribution.\n\n\n\n\n\n","category":"method"},{"location":"#StatisticalTests.lchoose-Tuple{Any,Any}","page":"Home","title":"StatisticalTests.lchoose","text":"lchoose(n, j)\n\nReturn the log of the number of choose j from n\n\n\n\n\n\n","category":"method"},{"location":"#StatisticalTests.t_test","page":"Home","title":"StatisticalTests.t_test","text":"t_test(X, y, paired::Bool=false, varequal::Bool=false)\n\nAn alias for fit(TTestModel, X, y, levels)\n\nThe arguments X and y can be a Matrix and a Vector or a Formula and DataFrame.\n\nThe argument type must be either :one, :paired, :simple, or welch.\n\nThe argument levels can be a Vector{String}.\n\nThe keyword argument μ0 can be a Real to use one sample t-test (default = 0.0). Ignored if the type is not :one.\n\n\n\n\n\n","category":"function"},{"location":"#StatisticalTests.t_test-Tuple{AbstractArray{T,1} where T<:Number,AbstractArray{T,1} where T<:Number}","page":"Home","title":"StatisticalTests.t_test","text":"t_test(x1, x2, paired::Bool=false, varequal::Bool=false)\n\nTwo sample t-test\n\nThe arguments x1 and x2 can be numeric Vector.\n\nPaired t-test will be performed if the keyword argument paired is true. Simple two sample t-test will be performed if the keyword argument varequal is true. Otherwise, Welch two sample t-test will be performed.\n\n\n\n\n\n","category":"method"},{"location":"#StatisticalTests.t_test-Tuple{AbstractArray{T,1} where T<:Number}","page":"Home","title":"StatisticalTests.t_test","text":"t_test(x, μ0::AbstractFloat=0.0)\n\nOne sample t-test\n\nThe argument x can be numeric Vector.\n\nOne sample t-test will be performed to compare between the mean of x1 and μ0.\n\n\n\n\n\n","category":"method"}]
}
