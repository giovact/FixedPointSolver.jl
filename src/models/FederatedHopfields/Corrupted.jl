struct Corrupted <: FPModel
    params::NamedTuple{(:β,:α,:r,:βT,:γ,:y,:ν),Tuple{FT,FT,FT,FT,FT,Int, Bool}}
    O::NamedVec #vector of orderparameters -> M, Q
    Oconj::NamedVec #vector of conjugate orderparameters -> Qhat
    InnerCW::HomogeneousSpin1Half
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors
Corrupted() = Corrupted(.1, .1,1.0,.1,1.0,2,false,ones(4)) # dummy constructor for testing
Corrupted(β,α,r,βT,γ,y,ν) = Corrupted(β,α,r,βT,γ,y,ν,zeros(4))
Corrupted(β,α,r,βT,γ,y,ν,O) = Corrupted((β=β,α=α,r=r,βT=βT,γ=γ,y=y,ν=ν),NamedArray(O,[:M;:P;:Q;:T]),NamedArray(zeros(3),[:M;:Q;:T]), HomogeneousSpin1Half(y),Dict{String,FT}())

create_model(Model::Corrupted, params) = Corrupted(params,Model.O, Model.Oconj, Model.InnerCW,Dict{String,Float64}())

set_integrationmethod(Model::Corrupted,X::TI) where TI <: IntegrationMethod = TracedSpin1HalfIntegral(Model,X)

########################################################## SELF-CONSISTENT EQUATIONS #########################################################
DeltaT(Model::Corrupted) =  1 - Model.params[:βT]
Delta(Model::Corrupted) =   1 - Model.params[:β] + Model.params[:β] * Model.O[:Q] + Model.params[:β] * Model.params[:r]^2 * (Model.O[:T] - Model.O[:Q]) 
Deltay(Model::Corrupted) =  1 - Model.params[:β] + Model.params[:β] * Model.O[:Q] - Model.params[:β] * Model.params[:r]^2 * (Model.O[:T] - Model.O[:Q]) * (Model.params[:y] - 1)

is_free_energy_definite(Model::Corrupted) = (Delta(Model)>0) && (Deltay(Model)>0)

function Deltas!(Model::Corrupted)
    Model.aux["ΔT"] = DeltaT(Model)
    Model.aux["Δ"] = clamp(Delta(Model), 1e-30, Inf)
    Model.aux["Δy"] = clamp(Deltay(Model), 1e-30, Inf)
end
    

#### Conjugate order params SP
function exprOconj!(Model::Corrupted)

    @extract Model : O Oconj params aux
    Deltas!(Model)

    β   = params[:β]
    βT  = params[:βT]
    r   = params[:r]
    α   = params[:α]
    y   = params[:y]

    M = O[:M]
    Q = O[:Q]
    T = O[:T]

    ΔT = aux["ΔT"]
    Δ  = aux["Δ"]
    Δy = aux["Δy"]

    ####################################################
    # 1) m-hat
    ####################################################
    Oconj[:M] = (α * βT * β * r^2 * M) / (ΔT * Δy)

    kappa = (βT * r^2 * M^2) / ΔT 
    ####################################################
    # 2) q-hat
    ####################################################

    # ----- 1 / Δ^2 contribution
    q_Δ2 = (α*β^2 / y^2) * (y - 1) * Q * (1-r^2)^2 / ( Δ^2 )

    # ----- 1 / Δy^2 contribution
    q_Δy2 = (α*β^2 / ( y^2 * Δy^2) ) * ( (1 + (y - 1) * r^2) * ( Q * (1 + r^2 * (y - 1)) + y * kappa )  ) 

    Oconj[:Q] = clamp( q_Δ2 + q_Δy2, 1e-10, Inf)

    ####################################################
    # 3) t-hat
    ####################################################

    # ----- contribution propto y (in r→ 1 gives original no noise eq)
    t_y = (β*r^2*y*β*r^2*(T-Q)) / (Δy*Δ) + β*r^2*y*( ((βT *β * r^2 * M^2) / ΔT ) + β*r^2*Q  )/(Δy^2)

    
    # ----- 1noise contribution (1- r^2)
    t_ΔyΔ2 = β*Q* β*r^2*(1-r^2) * ( (Δ^2-Δy^2)/(Δy^2*Δ^2)  )

    Oconj[:T] = clamp( (α / y) * (t_y + t_ΔyΔ2), Oconj[:Q] + 1e-10, Inf)


    # my re-writing
    
    # t_y = α * β^2 * r^4 * (T-Q) / (Δy*Δ)  +   α * β^2 * r^2 *( Q * r^2 + kappa )/( Δy^2) 
    # t_ΔyΔ2 = α * β^2 * r^2 * Q * (1 - r^2) * (Δ+Δy) * (β*r^2 * (T-Q) ) / ( Δ^2 * Δy^2 )
    # Oconj[:T] = clamp( t_y + t_ΔyΔ2, Oconj[:Q] + 1e-10, Inf)

    ##########################
    # check on auxiliary vars
    ##########################
    aux["Δ+Δy"] = Δ + Δy
    aux["Δ-Δy"] = Δ - Δy

end



function inner_trace(Model::Corrupted, z::FT)
    @extract Model : params O Oconj
    
    heff = Oconj[:M]  + sqrt(Oconj[:Q])*z
    Jeff = (params[:β]^params[:ν]) * params[:γ]/params[:y] + (Oconj[:T]-Oconj[:Q])

    return trace_CurieWeiss(Model.InnerCW,heff,Jeff)
end




function lhs!(Model::Corrupted,X::TracedSpin1HalfIntegral{TI},Onew::NamedVec) where TI <: IntegrationMethod

    exprOconj!(Model)
     
    @assert iszero(Onew)

    for s in [-1;1]
        inner_trace(Model, X)
        Onew[:M] += 0.5 * integrate(X.X,X.m)
        Onew[:P] += 0.0
        Onew[:Q] += 0.5 * integrate(X.X,X.m .^ 2)
        Onew[:T] += 0.5 * integrate(X.X,X.chi)
    end

    Onew[:P] = 0.0

    Onew[:Q] = clamp(Onew[:Q], 1 - 1 / Model.params[:β], 1.0)

    vareps = 1e-10
    Onew[:T] = clamp(Onew[:T], Onew[:Q] + vareps , 1.0)
    
end

########################################################## FREE ENERGY #########################################################

function free_energy(Model::Corrupted,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj aux 
    exprOconj!(Model)

    # energy_magnetizations = 0.5*O[:P]^2 + params[:α] * O[:M] * Oconj[:M]

    # energy_overlaps =  0.5*params[:α]*params[:β]*params[:y]*(O[:T]*Oconj[:T] - O[:Q]*Oconj[:Q] ) -  0.5*params[:α]*params[:β]*O[:T]*Oconj[:T]
    # logdetK =  (params[:α] / (2*params[:β])) * log(aux["Δ"]) +  (params[:α] / (2*params[:y]*params[:β])) * ( log(aux["Δy"]) - log(aux["Δ"]))
    # other = (0.5*params[:α] / aux["Δy"]) * (  O[:Q] + ( params[:βT] * O[:M]^2 / aux["ΔT"] ) )
    
    # ent =  0.0
    # for s in [-1;1]
    #     inner_trace(Model, X,s)
    #     ent += 0.5 * integrate(X.X,log.(X.Z))
    # end

    # return  energy_magnetizations + energy_overlaps + logdetK - other - ent/(params[:β] * params[:y])

    return -5
    

end


_isSG(Model::Corrupted,X::TI,Oold::NamedVec,tol::Float64) where TI <: IntegrationMethod = Model.O[:Q] < 1e-5 && Model.O[:P] < 1e-5
_isSignalR(Model::Corrupted,X::TI,Oold::NamedVec,tol::Float64) where TI <: IntegrationMethod = Model.O[:M] > 1e-2
_isExampleR(Model::Corrupted,X::TI,Oold::NamedVec,tol::Float64) where TI <: IntegrationMethod = Model.O[:P] > 1e-4 && Model.O[:M] < 1e-5
_isParamagnetic(Model::Corrupted,X::TI,Oold::NamedVec,tol::Float64) where TI <: IntegrationMethod = Model.O[:M] < 1e-3


#_is_RS_stable(Model::Corrupted,X::TI,Oold::NamedVec,tol::Float64) where TI <: IntegrationMethod  = min(Model.aux["replicon_parallel_parallel"], Model.aux["replicon_perpernd_perpernd"]) > 0 
_is_RS_stable(Model::Corrupted,X::TI,Oold::NamedVec,tol::Float64) where TI <: IntegrationMethod  = Model.aux["max_Lambda"] < 1.0

function eigs_dat!(Model::Corrupted,X::TracedSpin1HalfIntegral{TI}) where TI <: IntegrationMethod
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
    
    aux["replicon_parallel_parallel"] = replicon_parallel_parallel 
    aux["replicon_perpernd_perpernd"] = replicon_perpernd_perpernd
    aux["replicon_parallel_perpend"] = replicon_parallel_perpend

    aux["is_CS_ok"]  = aux["replicon_parallel_perpend"] > min(aux["replicon_parallel_parallel"], aux["replicon_perpernd_perpernd"])
    
    aux["is_RS_stable"] = min(aux["replicon_parallel_parallel"], aux["replicon_perpernd_perpernd"]) > 0 
    
end


function compute_aux!(Model::Corrupted,X::TracedSpin1HalfIntegral{TI}) where TI <: IntegrationMethod
    @extract Model : aux  params

#    eigs_dat!(Model, X)
    eigs_dat_NEW!(Model,X)
end




function eigs_dat_NEW!(Model::Corrupted,X::TracedSpin1HalfIntegral{TI}) where TI <: IntegrationMethod
    @extract Model : aux  params
    exprOconj!(Model)     
    inner_trace(Model, X)


    prefactor = (params[:α] * params[:β]^2) / params[:y]^2 
    aux["c_∥"] = 1 + (params[:y]-1)*params[:r]^2
    aux["c_⟂"] = 1-params[:r]^2
    aux["g"] = params[:y] -1 + params[:r]^2

    aux["k_∥∥"] = integrate(X.X, (1 .+ (params[:y]-1).*X.chi .- params[:y]*(X.m.^2) ).^2 )    
    aux["k_⟂⟂"] = integrate(X.X, (1 .- X.chi).^2 )
    aux["k_∥⟂"] = integrate(X.X, (1 .- X.chi) .* (1 .+ (params[:y]-1).*X.chi .- params[:y]*(X.m.^2) ) )
    
    # symmetric sector
    aux["Ts_11"] = prefactor * aux["k_∥∥"] * ( (aux["c_∥"] / aux["Δy"])^2 + (params[:y]-1)*(aux["c_⟂"] / aux["Δ"])^2)
    aux["Ts_12"] = prefactor * aux["k_∥∥"] * (params[:y]-1)*aux["c_⟂"] * ( (aux["c_∥"] / (aux["Δy"]^2)) + (aux["g"] / (aux["Δ"]^2)) )
    aux["Ts_21"] = prefactor * aux["k_⟂⟂"] * aux["c_⟂"] * ( (aux["c_∥"] / (aux["Δy"]^2)) + (aux["g"] / (aux["Δ"]^2)) )
    aux["Ts_22"] = prefactor * aux["k_⟂⟂"] * ( (params[:y]-1)*(aux["c_⟂"] / aux["Δy"])^2 + (aux["g"] / aux["Δ"])^2)
     

    eigs_Lambda_S = eigen([aux["Ts_11"] aux["Ts_12"]; aux["Ts_21"] aux["Ts_22"]]).values
    
    @assert isreal(eigs_Lambda_S[1])
    @assert isreal(eigs_Lambda_S[2])

    aux["Lambda_S"] = maximum(eigs_Lambda_S)

    

    ## A -sector
    
    aux["a11"] = (2 +(params[:y]-2)*params[:r]^2) / params[:y]
    aux["a12"] = (params[:y]-2)*aux["c_⟂"] / params[:y]^2
    aux["a21"] = 2*aux["c_⟂"] 
    aux["a22"] = (params[:y]-2+2*params[:r]^2) / params[:y]
    
    aux["e1_inv"] = aux["Δ"]*aux["Δy"]
    aux["e2_inv"] = aux["Δ"]^2
    
    aux["Ta_11"] = params[:α]*params[:β]^2 * aux["k_∥⟂"] * ( (aux["a11"]^2 / aux["e1_inv"]) + (aux["a12"]*aux["a21"]/aux["e2_inv"]))
    aux["Ta_12"] = params[:α]*params[:β]^2 * aux["k_∥⟂"] * aux["a12"] * ( (aux["a11"] / aux["e1_inv"]) + (aux["a22"]/aux["e2_inv"]))
    aux["Ta_21"] = params[:α]*params[:β]^2 * aux["k_⟂⟂"] * aux["a21"] * ( (aux["a11"] / aux["e1_inv"]) + (aux["a22"]/aux["e2_inv"]))
    aux["Ta_22"] = params[:α]*params[:β]^2 * aux["k_⟂⟂"] * ( (aux["a12"]*aux["a21"]/aux["e1_inv"]) + (aux["a22"]^2 / aux["e2_inv"]) )
    
    eigs_Lambda_A = eigen([aux["Ta_11"] aux["Ta_12"]; aux["Ta_21"] aux["Ta_22"]]).values
    
    @assert isreal(eigs_Lambda_A[1])
    @assert isreal(eigs_Lambda_A[2])


    if params[:y]==2
        aux["Lambda_A"] = params[:α] * params[:β]^2 *aux["k_∥⟂"] / (aux["Δ"]*aux["Δy"])
    else
        aux["Lambda_A"] = maximum(eigs_Lambda_A)
    end

    ## Y -sector
    if params[:y]>=4
        aux["Lambda_Y"] = params[:α] * params[:β]^2 * params[:r]^4 * aux["k_⟂⟂"] / aux["Δ"]^2
    else
        aux["Lambda_Y"] = -Inf
    end


    aux["max_Lambda"] = max(aux["Lambda_S"], aux["Lambda_A"], aux["Lambda_Y"])
    aux["is_RS_stable"] =  aux["max_Lambda"] < 1.0 
    
end


