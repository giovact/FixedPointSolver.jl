struct FederatedHopfield2SameAlpha <: FPModel
    params::NamedTuple{(:β,:α,:γ,:βT,:ν),Tuple{FT,FT,FT,FT,Int}}
    O::NamedVec #vector of orderparameters -> M
    Oconj::NamedVec #vector of orderparameters -> M
    aux::Dict{String,FT} # other stuff, among them free energy
end

#
FederatedHopfield2SameAlpha() = FederatedHopfield2SameAlpha(.8,0.1,.0,0.8,0,ones(4)*0.3) # dummy constructor for testing
FederatedHopfield2SameAlpha(β,α,γ,βT,ν) = FederatedHopfield2SameAlpha(β,α,γ,βT,ν,ones(4)) 
FederatedHopfield2SameAlpha(β,α,γ,βT,ν,O) = FederatedHopfield2SameAlpha((β=β,α=α,γ=γ,βT=βT,ν=ν),NamedArray(O,[:M1,:M2,:Q1,:Q2]),NamedArray(zeros(4),[:M1,:M2,:Q1,:Q2]), Dict{String, FT}())

create_model(Model::FederatedHopfield2SameAlpha, params) = FederatedHopfield2SameAlpha(params, Model.O,Model.Oconj,Dict{String, FT}())

########################################################## SELF-CONSISTENT EQUATIONS #########################################################


function aux_variables!(Model::FederatedHopfield2SameAlpha)
    @extract Model : O Oconj params aux 
    aux["ΔT"] = 1 - params[:βT]
    aux["Δ1"] = 1 - params[:β] * (1-O[:Q1])
    aux["Δ2"] = 1 - params[:β] * (1-O[:Q2])
end


function exprOconjall!(Model::FederatedHopfield2SameAlpha)
    @extract Model : O Oconj params aux
    aux_variables!(Model)

    Oconj[:M1] = ( params[:α] * params[:βT] * params[:β] * O[:M1] )  / (aux["ΔT"] * aux["Δ1"] )
    Oconj[:M2] = ( params[:α] * params[:βT] * params[:β] * O[:M2] )  / (aux["ΔT"] * aux["Δ2"] )

    
    Oconj[:Q1] = ( params[:α] * params[:βT] * (params[:β]^2) * (O[:M1]^2) )  / (aux["ΔT"] * aux["Δ1"]^2 ) + params[:α] * params[:β]^2 * O[:Q1]  /  aux["Δ1"]^2

    Oconj[:Q2] = ( params[:α] * params[:βT] * (params[:β]^2) * (O[:M2]^2) )  / (aux["ΔT"] * aux["Δ2"]^2 ) + params[:α] * params[:β]^2 * O[:Q2]  /  aux["Δ2"]^2

    
end


# order params SP

function internal_mag1(Model::FederatedHopfield2SameAlpha,z1,z2)
    @extract Model : params Oconj aux
    return tanh( Oconj[:M1] + sqrt(Oconj[:Q1])*z1  + atanh(tanh(0.5 * params[:γ] * params[:β]^params[:ν]) * tanh(Oconj[:M2] + sqrt(Oconj[:Q2])*z2) ) )
end


function internal_mag2(Model::FederatedHopfield2SameAlpha,z1,z2)
    @extract Model : params Oconj aux
    return tanh( Oconj[:M2] + sqrt(Oconj[:Q2])*z2  + atanh(tanh(0.5 * params[:γ] * params[:β]^params[:ν]) * tanh(Oconj[:M1] + sqrt(Oconj[:Q1])*z1) ) )
end

function exprM1(Model::FederatedHopfield2SameAlpha,X::TI) where TI <: IntegrationMethod
    IM1(z1,z2) = internal_mag1(Model,z1,z2)
    return integrate_double(X,IM1)
end

function exprM2(Model::FederatedHopfield2SameAlpha,X::TI) where TI <: IntegrationMethod
    @extract Model : params Oconj aux
    IM2(z1,z2) = internal_mag2(Model,z1,z2)
    return integrate_double(X,IM2)
end



function exprQ1(Model::FederatedHopfield2SameAlpha,X::TI) where TI <: IntegrationMethod
    IQ1(z1,z2) = internal_mag1(Model,z1,z2)^2
    return integrate_double(X,IQ1)
end

function exprQ2(Model::FederatedHopfield2SameAlpha,X::TI) where TI <: IntegrationMethod
    @extract Model : params Oconj aux
    IQ2(z1,z2) = internal_mag2(Model,z1,z2)^2
    return integrate_double(X,IQ2)
end



function lhs!(Model::FederatedHopfield2SameAlpha,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O
    exprOconjall!(Model)
    Onew[:M1] = exprM1(Model,X)
    Onew[:M2] = exprM2(Model,X)
    Onew[:Q1] = exprQ1(Model,X)
    Onew[:Q2] = exprQ2(Model,X)
end



########################################################## FREE ENERGY #########################################################
function free_energy(Model::FederatedHopfield2SameAlpha,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj aux 
    exprOconjall!(Model)

    piece_overlap = (0.5/params[:β])*Oconj[:Q1] * (1-O[:Q1]) + (0.5/params[:β])*Oconj[:Q2] * (1-O[:Q2])
    piece_magnetizations = O[:M1] * Oconj[:M1] / params[:β] + O[:M2] * Oconj[:M2] / params[:β]
    log_det = 0.5 * (params[:α] / params[:β]) * (log(aux["ΔT"]) + log(aux["Δ1"])) + 0.5 * (params[:α] / params[:β]) * (log(aux["ΔT"]) + log(aux["Δ2"])) 

    other_1 = 0.5 * params[:α] * O[:Q1] / aux["Δ1"] + 0.5 * params[:α] * (params[:βT] * O[:M1]^2) / (aux["ΔT"] * aux["Δ1"])

    other_2 = 0.5 * params[:α] * O[:Q2] / aux["Δ2"] + 0.5 * params[:α] * (params[:βT] * O[:M2]^2) / (aux["ΔT"] * aux["Δ2"])

    entropy = 

    Ientropy(z1,z2) = log(sum( exp( ξ1*( Oconj[:M1] + sqrt(Oconj[:Q1])*z1 ) + ξ2*( Oconj[:M2] + sqrt(Oconj[:Q2])*z2 )  + 0.5 * params[:γ] * params[:β]^params[:ν]* ξ1 * ξ2) for ξ1 in [-1;1] for ξ2 in [-1;1]))

    entropy = integrate_double(X,Ientropy)
    return piece_overlap + piece_magnetizations + log_det - other_1 - other_2 - entropy/params[:β]

end


_isSignal(Model::FederatedHopfield2SameAlpha, X::TI, Oold::NamedVec,tol::FT) where TI <: IntegrationMethod = (Model.O[:M1] > 1e-1) & (Model.O[:M2] > 1e-1)

########################################################## RS (dAT) STABILITY #########################################################
# Same replicon analysis as FederatedHopfield2.jl at α1 = α2 = α (derivation: dAT_SDN/dAT_SDD.pdf). The helpers
# dat_sdd2! and twospin_chis are defined in FederatedHopfield2.jl, which models_list.jl includes before this file.
# The order parameters M1, M2, Q1, Q2 are kept separate, so the 2x2 operator also covers student-asymmetric solutions;
# on the symmetric one (M1 = M2, Q1 = Q2) its eigenvalues are Λ± = αβ²(A ± B)/Δ², and the dAT line is αβ²(A + B) = Δ².
# aux keys: see dat_sdd2!.
function eigs_dat!(Model::FederatedHopfield2SameAlpha, X::TI) where TI <: IntegrationMethod
    @extract Model : params Oconj aux
    exprOconjall!(Model)                                    # Oconj and Δ1, Δ2 at the current O
    β = params[:β]
    J = 0.5 * params[:γ] * β^params[:ν]
    dat_sdd2!(aux, X, β, params[:α], params[:α], J, Oconj, aux["Δ1"], aux["Δ2"])
end

# called by solve_SPequations! at convergence (after free_energy!)
compute_aux!(Model::FederatedHopfield2SameAlpha, X::TI) where TI <: IntegrationMethod = eigs_dat!(Model, X)

# keepgoing function for critical_line.jl (--function is_RS_stable): stops when RS becomes unstable (dAT line)
_is_RS_stable(Model::FederatedHopfield2SameAlpha, X::TI, Oold::NamedVec, tol::FT) where TI <: IntegrationMethod = Model.aux["Lambda_dAT"] < 1.0

# P -> SG line at α1 = α2 = α: (T - 1)² = α (1 + tanh²J), J = γ β^ν / 2 (T_PSG in FederatedHopfield2.jl)
T_PSG(Model::FederatedHopfield2SameAlpha) = T_PSG(Model.params[:α], Model.params[:α], Model.params[:γ], Model.params[:ν])
