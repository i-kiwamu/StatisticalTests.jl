mutable struct AnovaResult<:TestResult
    type::Symbol
    contrast::Dict
    levels::Vector{String}
    DF::IntVector
    SS::FloatVector
    MSS::FloatVector
    Fval::FloatVector
    Pval::FloatVector
    ω²::FloatVector
end

function show(io::IO, pttr::AnovaResult)
end

function coeftable(pttr::AnovaResult)
end

"""
"""
function anova(
    formula::FormulaTerm{Term,Term},
    df::DataFrame,
    contrasts::Union{Dict,StatsModels.AbstractContrasts}=EffectsCoding(),
    method::Symbol=:auto,  # :auto, :ML, or :REML
    type::Symbol=:typeIII  # :typeI, :typeII, or :typeIII
)::AnovaResult
    # model type (LinearModel or MixedModel)
    if method == :REML
        model_type = MixedModel
    elseif method == :ML
        #
    elseif method == :auto
        #
    else
        error("Only :auto, :ML, or :REML are allowed for method argument.")
    end

    # specify contrasts & fit
    if contrasts isa Dict
        mod = fit(model_type, formula, df, contrasts=contrasts)
    else  # use the same contrast for all the categorical factors
        f = apply_schema(fomula, schema(formula, df))
        dict_contr = Dict()
        for term in f.rhs.terms
            if term isa StatsModels.CategoricalTerm
                get!(dict_contr, term.sym, contrasts)
            end
        end
        model_frame = ModelFrame(formula, df, contrasts=dict_contr)
        X = ModelMatrix(model_frame).m
        y = response(model_frame)
        mod = fit(model_type, X, y)
    end

    # specify type (I, II, III) for the sum of squares

end
