struct HopfieldTS <: FPModel
    params::NamedTuple{(:β,:α,:r,:βT),Tuple{FT,FT,FT,FT}}
    O::NamedVec #vector of orderparameters -> M, Q
    Oconj::NamedVec #vector of conjugate orderparameters -> Qhat
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors
HopfieldTS() = HopfieldTS(.1, .1,1.0,.1,ones(3)) # dummy constructor for testing
HopfieldTS(β,α,r,βT) = HopfieldTS(β,α,r,βT,zeros(3))
HopfieldTS(β,α,r,βT,O) = HopfieldTS((β=β,α=α,r=r,βT=βT),NamedArray(O,[:M;:P;:Q]),NamedArray(zeros(2),[:M;:Q]), Dict{String,FT}())

create_model(Model::HopfieldTS, params) = HopfieldTS(params,Model.O, Model.Oconj, Dict{String,FT}())


function sanity_checks!(Model::HopfieldTS)
    Model.params[:βT] >1 && throw("You should use a βT < 1 !!!!")
end

########################################################## SELF-CONSISTENT EQUATIONS #########################################################
DeltaS(Model::HopfieldTS) =  1 - Model.params[:β] * (1 - Model.O[:Q])
DeltaT(Model::HopfieldTS) =  1 - Model.params[:βT]
function Deltas!(Model::HopfieldTS)   
    Model.aux["ΔS"] = clamp(DeltaS(Model), 1e-30, Inf)
    Model.aux["ΔT"] = DeltaT(Model)
    Model.aux["Δ"] = Model.aux["ΔS"] * Model.aux["ΔT"]
end
    
is_free_energy_definite(Model::HopfieldTS) = Model.aux["ΔS"] >0 && Model.aux["ΔT"] >0

#### Conjugate order params SP
function exprOconj!(Model::HopfieldTS)
    @extract Model : O Oconj params aux
    Deltas!(Model)
    #Qhat
    Oconj[:M] = params[:βT] * params[:r]^2 * O[:M] / aux["Δ"]
    Oconj[:Q] = (params[:βT] * params[:r]^2 * O[:M]^2 + aux["ΔT"] * O[:Q] ) / (aux["ΔT"] * aux["ΔS"]^2)
end

χI(Model::HopfieldTS,z,s) = Model.params[:β] * ( Model.O[:P] * s + Model.params[:α] * Model.Oconj[:M] + sqrt(Model.params[:α] * Model.Oconj[:Q]) * z )

function exprM(Model::HopfieldTS,X::TI,s::Int) where TI <: IntegrationMethod
    IM(z) = tanh( χI(Model,z,s) )
    integrate(X,IM)
end

function exprQ(Model::HopfieldTS,X::TI,s::Int) where TI <: IntegrationMethod
    IQ(z) = tanh( χI(Model,z,s) )^2
    integrate(X,IQ)
end


function lhs!(Model::HopfieldTS,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O Oconj 
    exprOconj!(Model)
    @assert iszero(Onew)

    for s in [-1;1]
        m = exprM(Model, X, s)
        Onew[:M] += 0.5 * m
        Onew[:P] += 0.5 * s * m
        Onew[:Q] += 0.5 * exprQ(Model,X,s)
    end
    Onew[:P] = 0.0

end

########################################################## FREE ENERGY #########################################################
function exprEntropy(Model::HopfieldTS,X::TI,s::Int) where TI <: IntegrationMethod
    IQ(z) = log2cosh( χI(Model,z,s) )
    integrate(X,IQ)
end

function free_energy(Model::HopfieldTS,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj aux 
    exprOconj!(Model)

    energy_overlaps =  0.5*params[:α]*params[:β]*(1-O[:Q])*Oconj[:Q] 
    energy_magnetizations = 0.5 * O[:P]^2 + params[:α] * O[:M] * Oconj[:M]

    logdetK = (params[:α]/(2*params[:β])) * ( log(aux["ΔS"]) +  log(aux["ΔT"]) )

    other =  0.5*params[:α]* ( (O[:Q] / aux["ΔS"]) + (params[:βT]*params[:r]^2*O[:M]^2 / aux["Δ"]) )  

    entropy = 0.5 * ( exprEntropy(Model,X,1) + exprEntropy(Model,X,-1) )
    
    return  energy_magnetizations + energy_overlaps + logdetK - other -  entropy/params[:β]

end

_isSignalR(Model::HopfieldTS,X::TI,Oold::NamedVec,tol::Float64) where TI <: IntegrationMethod = Model.O[:M] > 1e-3
_isSignalR_1em3(Model::HopfieldTS,X::TI,Oold::NamedVec,tol::Float64) where TI <: IntegrationMethod = Model.O[:M] > 1e-3


function expr_dAT(Model::HopfieldTS,X::TI,j::Int) where TI <: IntegrationMethod
    I_dAT(z) = 1 / cosh( χI(Model,z,j) )^4
    integrate(X,I_dAT)
end


function eig_DAT(Model::HopfieldTS, X::TI) where TI <: IntegrationMethod
    @extract Model: O Oconj params aux
    exprOconj!(Model)
    λ_qblock= 0.0
    for s in [-1;1]
        λ_qblock += 0.5 * expr_dAT(Model,X,s)
    end
    λ_qblock = expr_dAT(Model,X,1)
    return aux["ΔS"]^2 - ( (params[:α]*params[:β]^2) ) *λ_qblock
end
_is_RS_stable(Model::HopfieldTS, X::TI, Oold::NamedVec,tol::FT) where TI <: IntegrationMethod = eig_DAT(Model,X)>0


function compute_aux!(Model::HopfieldTS,X::TI) where TI <: IntegrationMethod
    Model.aux["eig_DAT"] = eig_DAT(Model,X)
    Model.aux["is_RS_stable"] = Model.aux["eig_DAT"] > 0
    # instability of an SG fixed point (M = 0) towards the signal: loop gain λm = α ϱ X/(1-X) with X = β(1-Q) the SG
    # susceptibility and ϱ = βT r²/(1-βT); λm > 1 means that M = 0 is unstable (only meaningful at M = 0)
    Model.aux["X"] = 1 - DeltaS(Model)
    Model.aux["λm"] = Model.params[:α] * Model.params[:βT] * Model.params[:r]^2 / DeltaT(Model) * Model.aux["X"] / (1 - Model.aux["X"])
    Model.aux["is_SG_unstable"] = Model.aux["λm"] > 1
    Model.aux
end

# true while the fixed point is a stable SG (keepgoing function for critical_line.jl: scanning α at fixed β stops at the SG -> signal line)
_isSGstable(Model::HopfieldTS,X::TI,Oold::NamedVec,tol::FT) where TI <: IntegrationMethod = Model.O[:M] < 1e-3 && Model.aux["λm"] < 1

########################################################## SG -> SIGNAL LINE (RS, one scalar equation) #########################################################
# Below the triple point the SG (M = 0) becomes unstable towards the signal when α ϱ X/(1-X) = 1, with X = β(1-Q) the SG
# susceptibility and ϱ = βT r²/(1-βT). Imposing X = 1/(1+αϱ) fixes Q = 1 - X T and ΔS = 1 - X, so the RS SG equation becomes a
# single equation in (α, T). It is written as X = β E_z sech²(β h z) = (1/h) ∫ du φ(u T/h) sech²(u), h = sqrt(αQ)/ΔS, which stays
# smooth down to T = 0. Limits: α = 1/ϱ² at T_c = 1 + 1/ϱ (triple point), α = π/(2ϱ²) at T = 0, α ≈ α_c [1 + (2/3)(1 - T/T_c)²]
# near T_c. Derivation: dAT_SDN/HopfieldTS_1RSB/HopfieldTS_1RSB.pdf.
function residual_SGsignal_RS(α, T, ϱ)
    X = 1 / (1 + α * ϱ)
    h = sqrt(α * (1 - X * T)) / (1 - X)
    I = T == 0 ? 2 / sqrt(2π) : quadgk(u -> exp(-0.5 * (u * T / h)^2) / (sqrt(2π) * cosh(u)^2), -Inf, Inf)[1]
    return I / h - X
end

# α of the SG -> signal line at temperature T (NaN above the triple point, where there is no SG phase)
function alpha_SGsignal_RS(T, r, βT)
    ϱ = βT * r^2 / (1 - βT)
    T >= 1 + 1/ϱ && return NaN
    lo, hi = (1 + 1e-12) / ϱ^2, 2 / ϱ^2
    while residual_SGsignal_RS(hi, T, ϱ) * residual_SGsignal_RS(lo, T, ϱ) > 0
        hi *= 2
    end
    return find_zero(α -> residual_SGsignal_RS(α, T, ϱ), (lo, hi), Roots.Brent(); xatol = 1e-12)
end
alpha_SGsignal_RS(Model::HopfieldTS) = alpha_SGsignal_RS(1 / Model.params[:β], Model.params[:r], Model.params[:βT])

