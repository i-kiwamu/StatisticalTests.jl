mutable struct PairwiseTTestResult<:TestResult
    levels::Vector{String}
    correction::Symbol
    n_obs::IntVector
    means::FloatVector
    sds::FloatVector
    statistic::FloatVector
    dof::FloatVector
    pval::FloatVector
    cohen_d::FloatVector
end

function show(io::IO, pttr::PairwiseTTestResult)
end

function coeftable(pttr::PairwiseTTestResult)
end
