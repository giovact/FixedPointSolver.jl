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

_isSignalR(Model::HopfieldTS,X::TI,Oold::NamedVec,tol::Float64) where TI <: IntegrationMethod = Model.O[:M] > 1e-2
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
end

