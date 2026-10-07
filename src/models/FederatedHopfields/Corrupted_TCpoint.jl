struct Corrupted_TCpoint <: FPModel
    params::NamedTuple{(:r,:βT,:γ,:y,:ν),Tuple{FT,FT,FT,Int, Bool}}
    O::NamedVec #vector of orderparameters -> M, P,Q,T
    Oconj::NamedVec #vector of conjugate orderparameters -> Qhat
    InnerCW::HomogeneousSpin1Half
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors
Corrupted_TCpoint() = Corrupted_TCpoint(1.0,0.8,1.0, 2,false,ones(2)) # dummy constructor for testing
Corrupted_TCpoint(r,βT,γ,y,ν) = Corrupted_TCpoint(r,βT,γ,y,ν, ones(2))
Corrupted_TCpoint(r,βT,γ,y,ν,O) = Corrupted_TCpoint((r=r,βT=βT,γ=γ,y=y,ν=ν),NamedArray(O,[:βc,:αc]),NamedArray([1.0;0.0],[:T,:Th]), HomogeneousSpin1Half(y), Dict{String, FT}())
########################################################## SELF-CONSISTENT EQUATIONS #########################################################

create_model(Model::Corrupted_TCpoint, params) = Corrupted_TCpoint(params,Model.O, Model.Oconj, HomogeneousSpin1Half(Int(params[:y])),Dict{String,Float64}())

set_integrationmethod(Model::Corrupted_TCpoint,X::TI) where TI <: IntegrationMethod = TracedSpin1HalfIntegral(Model,X)

########################################################## SELF-CONSISTENT EQUATIONS #########################################################
DeltaT(Model::Corrupted_TCpoint) =  1 - Model.params[:βT]
Delta(Model::Corrupted_TCpoint) =  1 - Model.O[:βc] * ( 1 -  Model.params[:r]^2 * Model.Oconj[:T] )
Deltay(Model::Corrupted_TCpoint) =  1 - Model.O[:βc] * (1 + Model.params[:r]^2 *(Model.params[:y] - 1)*Model.Oconj[:T] ) 

is_free_energy_definite(Model::Corrupted_TCpoint) = (Delta(Model)>0) && (Deltay(Model)>0)

function Deltas!(Model::Corrupted_TCpoint)
    @extract Model : params O Oconj aux

    Model.aux["ΔT"] = DeltaT(Model)
    Model.aux["Δ"] = clamp(Delta(Model), 1e-30, Inf)
    Model.aux["Δy"] = clamp(Deltay(Model), 1e-30, Inf)

    Model.aux["Cq"] = ( (params[:y]-1) * (1-params[:r]^2)^2 / aux["Δ"]^2 ) + ( (1 + params[:r]^2 * (params[:y]-1) ) ^2 / aux["Δy"]^2 )
end

expr_chi(Model::Corrupted_TCpoint) = 1 + (Model.params[:y]-1)*Model.Oconj[:T]

function expr_OconjT!(Model::Corrupted_TCpoint) 
    @extract Model : params O Oconj aux
    Deltas!(Model)

    Oconj[:Th] = O[:αc] * O[:βc] * params[:r]^2 * Oconj[:T]  / ( aux["Δ"] * aux["Δy"] ) 
    Oconj[:T] = trace_CurieWeiss(Model.InnerCW,0.0,(Model.O[:βc]^params[:ν]) * params[:γ]/params[:y] + Oconj[:Th])[3]
    aux["chi"] = 1 + (params[:y]-1) * Oconj[:T]
end

function expr_betac(Model::Corrupted_TCpoint) 
    @extract Model : params O Oconj aux
    return params[:βT] * params[:r]^2 * params[:y]^2  / (aux["Cq"] * aux["ΔT"] * aux["Δy"] *  aux["chi"] )
end

function expr_alphac(Model::Corrupted_TCpoint) 
    @extract Model : params O Oconj aux
    return aux["ΔT"] * aux["Δy"] / ( params[:βT] * params[:r]^2 *  aux["chi"] * O[:βc] )
end

function lhs!(Model::Corrupted_TCpoint,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    expr_OconjT!(Model)
    Deltas!(Model)
    
    Onew[:βc] = expr_betac(Model)
    Onew[:αc] = expr_alphac(Model)

end

########################################################## FREE ENERGY #########################################################

function free_energy(Model::Corrupted_TCpoint,X::TI) where TI <: IntegrationMethod
    return -1.
end
