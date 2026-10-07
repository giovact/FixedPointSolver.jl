struct CoupledHopfieldTS <: FPModel
    params::NamedTuple{(:β,:α,:r,:βT,:γ,:y,:ν),Tuple{FT,FT,FT,FT,FT,Int, Bool}}
    O::NamedVec #vector of orderparameters -> M, Q
    Oconj::NamedVec #vector of conjugate orderparameters -> Qhat
    InnerCW::HomogeneousSpin1Half
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors
CoupledHopfieldTS() = CoupledHopfieldTS(.1, .1,1.0,.1,1.0,2,false,ones(4)) # dummy constructor for testing
CoupledHopfieldTS(β,α,r,βT,γ,y,ν) = CoupledHopfieldTS(β,α,r,βT,γ,y,ν,zeros(4))
CoupledHopfieldTS(β,α,r,βT,γ,y,ν,O) = CoupledHopfieldTS((β=β,α=α,r=r,βT=βT,γ=γ,y=y,ν=ν),NamedArray(O,[:M;:P;:Q;:T]),NamedArray(zeros(3),[:M;:Q;:T]), HomogeneousSpin1Half(y),Dict{String,FT}())

create_model(Model::CoupledHopfieldTS, params) = CoupledHopfieldTS(params,Model.O, Model.Oconj, Model.InnerCW,Dict{String,FT}())

set_integrationmethod(Model::CoupledHopfieldTS,X::TI) where TI <: IntegrationMethod = TracedSpin1HalfIntegral(Model,X)

########################################################## SELF-CONSISTENT EQUATIONS #########################################################
DeltaT(Model::CoupledHopfieldTS) =  1 - Model.params[:βT]
Delta(Model::CoupledHopfieldTS) =  1 - Model.params[:β] * (1 - Model.O[:T])
Deltay(Model::CoupledHopfieldTS) =  1 - Model.params[:β] * (1 - Model.O[:Q]) - (Model.params[:y]-1) * Model.params[:β] * (Model.O[:T] - Model.O[:Q])

is_free_energy_definite(Model::CoupledHopfieldTS) = (Delta(Model)>0) && (Deltay(Model)>0)

function Deltas!(Model::CoupledHopfieldTS)
    Model.aux["ΔT"] = DeltaT(Model)
    Model.aux["Δ"] = clamp(Delta(Model), 1e-20, Inf)
    Model.aux["Δy"] = clamp(Deltay(Model), 1e-20, Inf)
end
    

#### Conjugate order params SP
function exprOconj!(Model::CoupledHopfieldTS)
    @extract Model : O Oconj params aux
    Deltas!(Model)
    #Qhat
    Oconj[:M] = params[:βT] * params[:r]^2 * O[:M] / (aux["Δy"] * aux["ΔT"])
    Oconj[:Q] = (params[:βT] * params[:r]^2 * O[:M]^2 + aux["ΔT"] * O[:Q] ) / (aux["ΔT"] * aux["Δy"]^2)
    
    Oconj[:Q] = clamp(Oconj[:Q], 1e-10, Inf)
    Oconj[:T] = Oconj[:Q] + (O[:T] - O[:Q]) / (aux["Δ"] * aux["Δy"] )

    Oconj[:T] = clamp(Oconj[:T], 1e-10, Inf)

end


# function inner_trace(Model::CoupledHopfieldTS, z::FT,σ::Int)
#     @extract Model : params O Oconj
    
#     heff = params[:β] * ( O[:P] * σ + params[:α] * Oconj[:M]  + sqrt(params[:α]*Oconj[:Q])*z )
#     Jeff = (params[:β]^params[:ν]) * params[:γ]/params[:y] + params[:α]*params[:β]^2 * (Oconj[:T]-Oconj[:Q])

#     return trace_CurieWeiss(Model.InnerCW,heff,Jeff)
# end


function inner_trace(Model::CoupledHopfieldTS, z::FT)
    @extract Model : params O Oconj
    
    heff = params[:β] * ( params[:α] * Oconj[:M]  + sqrt(params[:α]*Oconj[:Q])*z )
    Jeff = (params[:β]^params[:ν]) * params[:γ]/params[:y] + params[:α]*params[:β]^2 * (Oconj[:T]-Oconj[:Q])

    return trace_CurieWeiss(Model.InnerCW,heff,Jeff)
end


# function inner_trace(Model::CoupledHopfieldTS, X,s)
#     @assert hasfield(typeof(Model),:InnerCW)
#     clean!(X)
#     for (iz, z) in enumerate(X.X.points)
#         Z, m, χ = inner_trace(Model,z,s)
#         X.Z[iz] = Z
#         X.m[iz] = m
#         X.chi[iz] = χ
#     end
# end


function lhs!(Model::CoupledHopfieldTS,X::TracedSpin1HalfIntegral{TI},Onew::NamedVec) where TI <: IntegrationMethod

    exprOconj!(Model)
     
    @assert iszero(Onew)

    # for s in [-1;1]
    #     inner_trace(Model, X,s)
    #     Onew[:M] += 0.5 * integrate(X.X,X.m)
    #     Onew[:P] += 0.5 * s * integrate(X.X,X.m)
    #     Onew[:Q] += 0.5 * integrate(X.X,X.m .^ 2)
    #     Onew[:T] += 0.5 * integrate(X.X,X.chi)
    # end

    inner_trace(Model, X)

    Onew[:M] = integrate(X.X,X.m)
    Onew[:Q] = integrate(X.X,X.m .^ 2)
    Onew[:T] = integrate(X.X,X.chi)

    # vareps = 1e-10
    # Onew[:Q] = clamp(Onew[:Q], 1 - 1 / Model.params[:β], 1.0-vareps)

    # Onew[:T] = clamp(Onew[:T], Onew[:Q] + vareps , 1.0-vareps)



end




########################################################## FREE ENERGY #########################################################

function free_energy(Model::CoupledHopfieldTS,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj aux 
    exprOconj!(Model)

    energy_magnetizations = 0.5*O[:P]^2 + params[:α] * O[:M] * Oconj[:M]

    energy_overlaps =  0.5*params[:α]*params[:β]*params[:y]*(O[:T]*Oconj[:T] - O[:Q]*Oconj[:Q] ) -  0.5*params[:α]*params[:β]*O[:T]*Oconj[:T]
    logdetK =  (params[:α] / (2*params[:β])) * log(aux["Δ"]) +  (params[:α] / (2*params[:y]*params[:β])) * ( log(aux["Δy"]) - log(aux["Δ"]))
    other = (0.5*params[:α] / aux["Δy"]) * (  O[:Q] + ( params[:βT] * O[:M]^2 / aux["ΔT"] ) )
    
    ent =  0.0
    # for s in [-1;1]
    #     inner_trace(Model, X,s)
    #     ent += 0.5 * integrate(X.X,log.(X.Z))
    # end

    return  energy_magnetizations + energy_overlaps + logdetK - other - ent/(params[:β] * params[:y])
    

end


_isSG(Model::CoupledHopfieldTS,X::TI,Oold::NamedVec,tol::FT) where TI <: IntegrationMethod = Model.O[:Q] > 1e-4 && Model.O[:P] < 1e-5
_isSignalR(Model::CoupledHopfieldTS,X::TI,Oold::NamedVec,tol::FT) where TI <: IntegrationMethod = Model.O[:M] > 1e-2
_isExampleR(Model::CoupledHopfieldTS,X::TI,Oold::NamedVec,tol::FT) where TI <: IntegrationMethod = Model.O[:P] > 1e-4 && Model.O[:M] < 1e-5

_isParamagnetic(Model::CoupledHopfieldTS,X::TI,Oold::NamedVec,tol::FT) where TI <: IntegrationMethod = Model.O[:M] < 1e-4 

function _is_RS_stable(Model::CoupledHopfieldTS,X::TI,Oold::NamedVec,tol::Float64) where TI <: IntegrationMethod 
    replicon_parallel_parallel, replicon_perpernd_perpernd, replicon_parallel_perpend = eigs_dat(Model, X)

    return  min(replicon_parallel_parallel,replicon_perpernd_perpernd) >0
end 


function eigs_dat(Model::CoupledHopfieldTS,X::TI) where TI <: IntegrationMethod
    @extract Model : aux  params
    exprOconj!(Model)     
    inner_trace(Model, X)
    chi_parallel_parallel = integrate(X.X, (1 .- X.m) .^2) 
    chi_perpernd_perpernd = integrate(X.X, (1 .- X.chi) .^2)
    chi_parallel_perpend =  integrate(X.X, (1 .- X.m) .* (1 .- X.chi) )

    pre_factor = params[:α] * params[:β]^2 / params[:y]^2
    replicon_parallel_parallel = 1 - pre_factor * ( (1 + (params[:y] - 1)*params[:r]^2 )^2 / aux["Δy"]^2 ) * chi_parallel_parallel
    replicon_perpernd_perpernd = 1 - pre_factor * ( (1 - params[:r]^2 )^2 / aux["Δ"]^2 ) * chi_perpernd_perpernd
    replicon_parallel_perpend = 1 - pre_factor * ( ((1 + (params[:y] - 1)*params[:r]^2 ) * (1-params[:r]^2)) / (aux["Δ"] * aux["Δy"]) ) * chi_parallel_perpend
    
    return replicon_parallel_parallel, replicon_perpernd_perpernd, replicon_parallel_perpend
end



function compute_aux!(Model::CoupledHopfieldTS,X::TI) where TI <: IntegrationMethod
    @extract Model : aux  params

    replicon_parallel_parallel, replicon_perpernd_perpernd, replicon_parallel_perpend = eigs_dat(Model, X)
    aux["replicon_parallel_parallel"] = replicon_parallel_parallel 
    aux["replicon_perpernd_perpernd"] = replicon_perpernd_perpernd
    aux["replicon_parallel_perpend"] = replicon_parallel_perpend


    aux["is_CS_ok"]  = aux["replicon_parallel_perpend"] > min(aux["replicon_parallel_parallel"], aux["replicon_perpernd_perpernd"])
    
    aux["is_RS_stable"] = min(aux["replicon_parallel_parallel"], aux["replicon_perpernd_perpernd"]) > 0 
    
end
