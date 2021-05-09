mutable struct AnovaResult<:TestResult
    type::Symbol
    contrast::Dict
    levels::Vector{String}
    dof::IntVector
    ss::FloatVector
    mss::FloatVector
    fval::FloatVector
    pval::FloatVector
    omega2::FloatVector
end

function show(io::IO, pttr::AnovaResult)
end

function coeftable(pttr::AnovaResult)
end

"""
"""
function anova(
    mod::StatsModels.RegressionModel,
    type::Symbol=:auto
)::AnovaResult
end
