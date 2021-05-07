using StatisticalTests
using Test, RDatasets, Statistics, StatsModels, Distributions

test_show(x) = show(IOBuffer(), x)

include("test_t.jl")
include("test_f.jl")
include("test_ks.jl")
include("test_wilcox.jl")
