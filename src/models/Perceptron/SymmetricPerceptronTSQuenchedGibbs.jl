struct SymmetricPerceptronTSQuenchedGibbs <: FPModel
    params::NamedTuple{(:β,:α,:k),Tuple{Float64,Float64,Float64}}
    O::NamedVec #vector of orderparameters -> M, Q
    Oconj::NamedVec #vector of conjugate orderparameters -> Qhat
    aux::Dict{String,Float64} # other stuff, among them free energy
end

# constructors
SymmetricPerceptronTSQuenchedGibbs() = SymmetricPerceptronTSQuenchedGibbs(1.0,1.0,1.0,ones(2)*1e-2 ) # dummy constructor for testing
SymmetricPerceptronTSQuenchedGibbs(β,α,k) = SymmetricPerceptronTSQuenchedGibbs(β,α,k,zeros(2))
SymmetricPerceptronTSQuenchedGibbs(β,α,k,O) = SymmetricPerceptronTSQuenchedGibbs((β=β,α=α,k=k),NamedArray(O,[:R;:Q]),NamedArray(zeros(2),[:R;:Q]), Dict{String, Float64}())

create_model(Model::SymmetricPerceptronTSQuenchedGibbs, params) = SymmetricPerceptronTSQuenchedGibbs(params, Model.O, Model.Oconj, Dict{String, Float64}())

function aux_variables!(Model::SymmetricPerceptronTSQuenchedGibbs)
    Model.aux["QmR2"] = Model.O[:Q] - Model.O[:R]^2
    Model.aux["1mQ"] = 1 - Model.O[:Q]
end

term_A_plus(Model::SymmetricPerceptronTSQuenchedGibbs, t) = (Model.params[:k] * sqrt(Model.O[:Q]) + (Model.O[:R]) * t) / sqrt(2 * Model.aux["QmR2"])
term_A_minus(Model::SymmetricPerceptronTSQuenchedGibbs, t) = (Model.params[:k] * sqrt(Model.O[:Q]) - (Model.O[:R]) * t) / sqrt(2 * Model.aux["QmR2"])

term_B_plus(Model::SymmetricPerceptronTSQuenchedGibbs, t) = (Model.params[:k] + sqrt(Model.O[:Q]) * t) / sqrt(2 * Model.aux["1mQ"])
term_B_minus(Model::SymmetricPerceptronTSQuenchedGibbs, t) = (Model.params[:k] - sqrt(Model.O[:Q]) * t) / sqrt(2 * Model.aux["1mQ"])

num_A_plus(Model::SymmetricPerceptronTSQuenchedGibbs, t) = (Model.params[:k] * sqrt(Model.O[:Q]) + (Model.O[:R]) * t)
num_A_minus(Model::SymmetricPerceptronTSQuenchedGibbs, t) = (Model.params[:k] * sqrt(Model.O[:Q]) - (Model.O[:R]) * t)

num_A_plus_R(Model::SymmetricPerceptronTSQuenchedGibbs, t) = Model.O[:R] * (Model.params[:k] * sqrt(Model.O[:Q]) + (Model.O[:R]) * t)
num_A_minus_R(Model::SymmetricPerceptronTSQuenchedGibbs, t) = Model.O[:R] * (Model.params[:k] * sqrt(Model.O[:Q]) - (Model.O[:R]) * t)

num_B_plus(Model::SymmetricPerceptronTSQuenchedGibbs, t) = (Model.params[:k] + sqrt(Model.O[:Q]) * t)
num_B_minus(Model::SymmetricPerceptronTSQuenchedGibbs, t) = (Model.params[:k] - sqrt(Model.O[:Q]) * t)

# Conjugate order params SP
function exprOconjall!(Model::SymmetricPerceptronTSQuenchedGibbs, X::TI) where TI <: IntegrationMethod
    @extract Model params O Oconj aux
    aux_variables!(Model)

    diff_A_q(t) = (1 / sqrt(pi)) * (
    exp(-((term_A_minus(Model, t))^2)) * (Model.params[:k] / (2 * sqrt(2) * sqrt(Model.O[:Q]) * sqrt(Model.aux["QmR2"])) - num_A_minus(Model, t) / (2 * sqrt(2) * (Model.aux["QmR2"])^(3/2))) + 
    exp(-((term_A_plus(Model, t))^2)) * (Model.params[:k] / (2 * sqrt(2) * sqrt(Model.O[:Q]) * sqrt(Model.aux["QmR2"])) - num_A_plus(Model, t) / (2 * sqrt(2) * (Model.aux["QmR2"])^(3/2))))

    diff_B_q(t) = (1 / sqrt(pi)) * (
    exp(-((term_B_minus(Model, t))^2)) * (-(t / (2 * sqrt(2) * sqrt(Model.aux["1mQ"]) * sqrt(Model.O[:Q]))) + 
    num_B_minus(Model, t) / (2 * sqrt(2) * (Model.aux["1mQ"])^(3/2))) + 
    exp(-((term_B_plus(Model, t))^2)) * (t / (2 * sqrt(2) * sqrt(Model.aux["1mQ"]) * sqrt(Model.O[:Q])) + 
    num_B_plus(Model, t) / (2 * sqrt(2) * (Model.aux["1mQ"])^(3/2))))

    diff_A_R(t) = (1 / sqrt(pi)) * (
    exp(-((term_A_minus(Model, t))^2)) * (-(t / (sqrt(2) * sqrt(Model.aux["QmR2"]))) + 
    num_A_minus_R(Model, t) / (sqrt(2) * (Model.aux["QmR2"])^(3/2))) + 
    exp(-((term_A_plus(Model, t))^2)) * (t / (sqrt(2) * sqrt(Model.aux["QmR2"])) + 
    num_A_plus_R(Model, t) / (sqrt(2) * (Model.aux["QmR2"])^(3/2))))
    
    erf_A(t) = 0.5 * (erf(term_A_plus(Model, t)) + erf(term_A_minus(Model, t)))
    
    erfc_A(t) = 0.5 * (1 - erf(term_A_plus(Model, t))) + 0.5 * (1 - erf(term_A_minus(Model, t)))
    
    erf_B(t) = 0.5 * (erf(term_B_plus(Model, t)) + erf(term_B_minus(Model, t)))
    
    erfc_B(t) = 0.5 * (1 - erf(term_B_plus(Model, t))) + 0.5 * (1 - erf(term_B_minus(Model, t)))
    
    part1_q(t) = erfc_B(t) + exp(-Model.params[:β]) * erf_B(t)
    part2_q(t) = exp(-Model.params[:β]) * erfc_B(t) + erf_B(t)
    
    f1_q(t) = - diff_A_q(t) * log(part1_q(t))
    
    f2_q(t) = erfc_A(t) * ((- diff_B_q(t) + exp(-Model.params[:β]) * diff_B_q(t)) / part1_q(t))
    
    f3_q(t) = diff_A_q(t) * log(part2_q(t))
    
    f4_q(t) = erf_A(t) * ((- exp(-Model.params[:β]) * diff_B_q(t) + diff_B_q(t)) / part2_q(t))
    
    integrand_q = integrate(X,f1_q) + integrate(X,f2_q) + integrate(X,f3_q) + integrate(X,f4_q)
    
    f1_R(t) = - diff_A_R(t) * log(part1_q(t))
    
    f2_R(t) = diff_A_R(t) * log(part2_q(t))
    
    integrand_R = integrate(X,f1_R) + integrate(X,f2_R)
    
    Oconj[:R]  = params[:α] * integrand_R
    eps = 1e-15
    Oconj[:Q] = clamp(- 2 * params[:α] * integrand_q, eps, Inf)
end

# Order params SP
function exprR(Model::SymmetricPerceptronTSQuenchedGibbs, X::TI) where TI <: IntegrationMethod
    @extract Model : params Oconj O 
    O[:R] == 0.0 && return 0.0
    
    IR(z) = tanh(sqrt(Oconj[:Q]) * z + Oconj[:R])
    return integrate(X, z -> IR(z))
end

function exprQ(Model::SymmetricPerceptronTSQuenchedGibbs, X::TI) where TI <: IntegrationMethod
    @extract Model : params Oconj 
    
    IQ(z) = tanh(sqrt(Oconj[:Q]) * z + Oconj[:R]) ^ 2
    return integrate(X, z -> IQ(z))
end

function lhs!(Model::SymmetricPerceptronTSQuenchedGibbs, X::TI, Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model :O Oconj
    exprOconjall!(Model, X)
    
    Onew[:R] = exprR(Model, X)
    Onew[:Q] = exprQ(Model,X) #clamp(real(exprQ(Model,X)), Onew[:R]^2 + eps, 1 - eps)  
end

########################################################## FREE ENERGY #########################################################


function free_energy(Model::SymmetricPerceptronTSQuenchedGibbs,X::TI) where TI <: IntegrationMethod
    @extract Model params O Oconj aux
    exprOconjall!(Model,X)

    diff_A_q(t) = (1 / sqrt(pi)) * (
    exp(-((term_A_minus(Model, t))^2)) * (Model.params[:k] / (2 * sqrt(2) * sqrt(Model.O[:Q]) *                         sqrt(Model.aux["QmR2"])) - num_A_minus(Model, t) / (2 * sqrt(2) * (Model.aux["QmR2"])^(3/2))) + 
    exp(-((term_A_plus(Model, t))^2)) * (Model.params[:k] / (2 * sqrt(2) * sqrt(Model.O[:Q]) *                         sqrt(Model.aux["QmR2"])) - num_A_plus(Model, t) / (2 * sqrt(2) * (Model.aux["QmR2"])^(3/2))))

    diff_B_q(t) = (1 / sqrt(pi)) * (
    exp(-((term_B_minus(Model, t))^2)) * (-(t / (2 * sqrt(2) * sqrt(Model.aux["1mQ"]) * sqrt(Model.O[:Q]))) + 
    num_B_minus(Model, t) / (2 * sqrt(2) * (Model.aux["1mQ"])^(3/2))) + 
    exp(-((term_B_plus(Model, t))^2)) * (t / (2 * sqrt(2) * sqrt(Model.aux["1mQ"]) * sqrt(Model.O[:Q])) + 
    num_B_plus(Model, t) / (2 * sqrt(2) * (Model.aux["1mQ"])^(3/2))))

    diff_A_R(t) = (1 / sqrt(pi)) * (
    exp(-((term_A_minus(Model, t))^2)) * (-(t / (sqrt(2) * sqrt(Model.aux["QmR2"]))) + 
    num_A_minus_R(Model, t) / (sqrt(2) * (Model.aux["QmR2"])^(3/2))) + 
    exp(-((term_A_plus(Model, t)) ^2)) * (t / (sqrt(2) * sqrt(Model.aux["QmR2"])) + 
    num_A_plus_R(Model, t) / (sqrt(2) * (Model.aux["QmR2"])^(3/2))))
    
    erf_A(t) = 0.5 * (erf(term_A_plus(Model, t)) + erf(term_A_minus(Model, t)))
    
    erfc_A(t) = 0.5 * (1 - erf(term_A_plus(Model, t))) + 0.5 * (1 - erf(term_A_minus(Model, t)))
    
    erf_B(t) = 0.5 * (erf(term_B_plus(Model, t)) + erf(term_B_minus(Model, t)))
    
    erfc_B(t) = 0.5 * (1 - erf(term_B_plus(Model, t))) + 0.5 * (1 - erf(term_B_minus(Model, t)))
    
    part1_q(t) = erfc_B(t) + exp(-Model.params[:β]) * erf_B(t)
    part2_q(t) = exp(-Model.params[:β]) * erfc_B(t) + erf_B(t)
    
    f1_q(t) = - diff_A_q(t) * log(part1_q(t))
    
    f2_q(t) = erfc_A(t) * (((- 1 + exp(-Model.params[:β])) * diff_B_q(t)) / part1_q(t))
    
    f3_q(t) = diff_A_q(t) * log(part2_q(t))
    
    f4_q(t) = erf_A(t) * (((1 - exp(-Model.params[:β])) * diff_B_q(t)) / part2_q(t))
    
    integrand_q = integrate(X,f1_q) + integrate(X,f2_q) + integrate(X,f3_q) + integrate(X,f4_q)
    
    f1_R(t) = - diff_A_R(t) * log(part1_q(t))
    
    f2_R(t) = diff_A_R(t) * log(part2_q(t))
    
    integrand_R = integrate(X,f1_R) + integrate(X,f2_R)
    
    Oconj[:R]  = params[:α] * integrand_R
    eps = 1e-15
    Oconj[:Q] = clamp(- 2 * params[:α] * integrand_q, eps, Inf)
    

    constraints = 0.5 * Oconj[:Q] * (O[:Q] - 1) - O[:R] * Oconj[:R]
    
    I_GS(t) = log(2 * cosh(sqrt(Oconj[:Q])*t + Oconj[:R] ))
    
    int1(t) = erfc_A(t) * log(part1_q(t))
    int2(t) = erf_A(t) * log(part2_q(t))
    
    GE = integrate(X,int1) + integrate(X,int2)
    
    G = constraints + integrate(X,I_GS)  + params[:α] * GE
    
    return - G # this quantity is beta*f
end

αc_2ndorder_quenched(k,β) = 0.5 * exp(k^2) * π * (erf(sqrt(0.5)*k) + exp(-β)*erfc(sqrt(0.5)*k) ) * (erfc(sqrt(0.5)*k) + exp(-β)*erf(sqrt(0.5)*k) )  / ( k^2 * (1-exp(-2 * β))) 


_isR_notTeacher(Model::SymmetricPerceptronTSQuenchedGibbs,X::TI,Oold::NamedVec,tol::FT) where {TI <: IntegrationMethod} =  Model.O[:R] < 0.99
_isR_Teacher(Model::SymmetricPerceptronTSQuenchedGibbs,X::TI,Oold::NamedVec,tol::FT) where {TI <: IntegrationMethod} =  Model.O[:R] >= 0.99
_isf_negative(Model::SymmetricPerceptronTSQuenchedGibbs,X::TI,Oold::NamedVec,tol::FT) where {TI <: IntegrationMethod} =  (Model.aux["f"] < 0.0) & (_isR_notTeacher(Model, X, Oold, tol))
_isR_null(Model::SymmetricPerceptronTSQuenchedGibbs,X::TI,Oold::NamedVec,tol::FT) where {TI <: IntegrationMethod} =  Model.O[:R] < 1e-2
_isTeacher_dropping(Model::SymmetricPerceptronTSQuenchedGibbs,X::TI,Oold::NamedVec,tol::FT) where {TI <: IntegrationMethod} =  Oold[:R] - Model.O[:R] < 0.2


