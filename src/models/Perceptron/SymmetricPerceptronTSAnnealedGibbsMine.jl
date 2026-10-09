struct SymmetricPerceptronTSAnnealedGibbsMine <: FPModel # this is the structure name corresponding to our model. For coding reasons it is defined as a "subtype" of the type FPModel. So it's one subtype of a generic model you want to find the solution of a certain set of fixed point equations of
    params::NamedTuple{(:α,:k,:β),Tuple{FT,FT,FT}}
    O::NamedVec #vector of orderparameters -> in this case it will have just one element (R). Named Vectors can be accessed also by their names' components
    Oconj::NamedVec
    U::Interval
    RbutU::IntervalSet
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors
SymmetricPerceptronTSAnnealedGibbsMine() = SymmetricPerceptronTSAnnealedGibbsMine(1.0,1.0,1.0,ones(1)*1e-2) # dummy constructor for testing
SymmetricPerceptronTSAnnealedGibbsMine(α,k,β) = SymmetricPerceptronTSAnnealedGibbsMine(α,k,β,zeros(1)) # this constructs a model with parameters α, k, β and R initialized to 0

# main constructor 
SymmetricPerceptronTSAnnealedGibbsMine(α,k,β,O) = SymmetricPerceptronTSAnnealedGibbsMine((α=α,k=k,β=β), # this first line creates the named tuple corresponding to the field "params"
                                                                                NamedArray(O,[:R]), # this creates a named vector starting from the input vector O (which has 1 entry only) and it assigns the name of its unique component as "R"
                                                                                NamedArray(zeros(1),[:R]),
                                                                                Interval(-k,k),
                                                                                IntervalSet((Interval(-Inf,-k), Interval(k,Inf))),
                                                                                Dict{String, FT}() # empty dictionary initialization (then it will be filled with e.g. the free energy value)
                                                                                )

create_model(Model::SymmetricPerceptronTSAnnealedGibbsMine, params) = SymmetricPerceptronTSAnnealedGibbsMine(params, Model.O, Model.Oconj,Interval(-params[:k], params[:k]),  IntervalSet((Interval(-Inf,-params[:k]), Interval(params[:k],Inf))), Dict{String, FT}())


########################################################## SELF-CONSISTENT EQUATIONS #########################################################

function aux_variables!(Model::SymmetricPerceptronTSAnnealedGibbsMine)
    # if you need to compute quantities that you will use often in the calculations (e.g. 1-R^2), you can save them in the dictionary
    # for instance
    Model.aux["1mR2"] = 1 - Model.O[:R]^2 

end

#### Conjugate order params SP

# this is just an example of old functions I had used in the symmetric perceptron (with no teacher) to shorten the computation of all the arguments of the error functions  
Hp(Model::SymmetricPerceptronTSAnnealedGibbsMine,ω) = (Model.params[:k] + Model.O[:R] * ω ) / sqrt(Model.aux["1mR2"])
Hm(Model::SymmetricPerceptronTSAnnealedGibbsMine,ω) = (Model.params[:k] - Model.O[:R] * ω ) / sqrt(Model.aux["1mR2"])

Hp_dR(Model::SymmetricPerceptronTSAnnealedGibbsMine,ω) = (Model.params[:k] * Model.O[:R] +  ω ) / sqrt(( Model.aux["1mR2"]) ^3 )
Hm_dR(Model::SymmetricPerceptronTSAnnealedGibbsMine,ω) = (Model.params[:k] * Model.O[:R] -  ω ) / sqrt(( Model.aux["1mR2"]) ^3 )

Herf_der(x) = - exp(-0.5 * x^2) / sqrt(2π)
Herfbar(x) = 0.5 * erf(sqrt(0.5)*x)
Herfbar_der(x) = exp(-0.5 * x^2) / sqrt(2π)

function GE(Model::SymmetricPerceptronTSAnnealedGibbsMine,X::TI) where TI <: IntegrationMethod
     
    intGE_U(ω) = ω ∈ Model.U ?  Herfbar(Hp(Model,ω)) + Herfbar(Hm(Model,ω)) + exp(-Model.params[:β]) * ( Herf(Hp(Model,ω)) + Herf(Hm(Model,ω))  )  : 0.0

    intGE_RbutU(ω) = ω ∈ Model.RbutU ? exp(-Model.params[:β]) * ( Herfbar(Hp(Model,ω)) + Herfbar(Hm(Model,ω)) )  +  ( Herf(Hp(Model,ω)) + Herf(Hm(Model,ω))  ) : 0.0

    return integrate(X,intGE_U) + integrate(X,intGE_RbutU)
end


function exprRhat!(Model::SymmetricPerceptronTSAnnealedGibbsMine,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj aux
    aux_variables!(Model)

    intGEder_U(ω) =   ω ∈ Model.U ?  ( exp(- 0.5 * Hp(Model,ω)^2) * Hp_dR(Model, ω) / sqrt(2* π) ) +  ( exp(- 0.5 * Hm(Model,ω)^2) * Hm_dR(Model, ω) / sqrt(2* π) ) : 0.0
    intGEder_RbutU(ω) = ω ∈ Model.RbutU ?   ( exp(- 0.5 * Hp(Model,ω)^2) * Hp_dR(Model, ω) / sqrt(2* π) ) +  ( exp(- 0.5 * Hm(Model,ω)^2) * Hm_dR(Model, ω) / sqrt(2* π) ) : 0.0
    Oconj[:R] = params[:α] * (1 - exp(-params[:β])) * ( integrate(X,intGEder_U) - integrate(X,intGEder_RbutU)  ) / GE(Model, X)

end

# order params SP

function lhs!(Model::SymmetricPerceptronTSAnnealedGibbsMine,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O Oconj
    
    exprRhat!(Model,X)
    Onew[:R] = tanh(Oconj[:R])
end

########################################################## FREE ENERGY #########################################################

function free_energy(Model::SymmetricPerceptronTSAnnealedGibbsMine,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj
    exprRhat!(Model,X)
    
    G = - O[:R] * Oconj[:R] + params[:α] * log(GE(Model,X)) + log2cosh(Oconj[:R])


    return - G # this quantity is beta*f
end


function compute_aux!(Model::SymmetricPerceptronTSAnnealedGibbsMine,X::TI) where TI <: IntegrationMethod
    int_unif_U(ω) = ω ∈ Model.U ? 1.0 : 0.0

    int_unif_RbutU(ω) = ω ∈ Model.RbutU ? 1.0 : 0.0

    Model.aux["intU"] = integrate(X,int_unif_U)
    Model.aux["intRbutU"] = integrate(X,int_unif_RbutU)
end

αc_2ndorder(Model::SymmetricPerceptronTSAnnealedGibbsMine) = αc_2ndorder_annealed(Model.params[:k],Model.params[:β])
αc_2ndorder_annealed(k,β) = 0.25 * π  * exp(k^2)  * (erfc(sqrt(0.5)*k)^2 + erf(sqrt(0.5)*k)^2 + 2*erf(sqrt(0.5)*k) * erfc(sqrt(0.5)*k) * exp(-β)  ) / ( k^2 * (1-exp(-β))) 


αc_1storder(Model::SymmetricPerceptronTSAnnealedGibbsMine) = αc_1storder_annealed(Model.params[:k],Model.params[:β])
αc_1storder_annealed(k,β) = -log(2) / log(erfc(sqrt(0.5)*k)^2 + erf(sqrt(0.5)*k)^2 + 2*erf(sqrt(0.5)*k) * erfc(sqrt(0.5)*k) * exp(-β)  )

_isRnull(Model::SymmetricPerceptronTSAnnealedGibbsMine,X::TI,Oold::NamedVec,tol::FT) where {TI <: IntegrationMethod} =  Model.O[:R] < 1e-2

_isR_notTeacher(Model::SymmetricPerceptronTSAnnealedGibbsMine,X::TI,Oold::NamedVec,tol::FT) where {TI <: IntegrationMethod} =  Model.O[:R] < 0.99
_isR_notTeacher_stronger(Model::SymmetricPerceptronTSAnnealedGibbsMine,X::TI,Oold::NamedVec,tol::FT) where {TI <: IntegrationMethod} =  Model.O[:R] < 0.999


_isRsubopt_butdominant(Model::SymmetricPerceptronTSAnnealedGibbsMine,X::TI,Oold::NamedVec,tol::FT) where {TI <: IntegrationMethod} =  Model.O[:R] < 1e-3 || Model.aux["f"] < 0.0