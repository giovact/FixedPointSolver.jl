struct SymmetricPerceptronTSQuenchedGibbsfixedR <: FPModel
    params::NamedTuple{(:β,:α,:k,:R),Tuple{Float64,Float64,Float64,Float64}}
    O::NamedVec #vector of orderparameters -> M, Q
    Oconj::NamedVec #vector of conjugate orderparameters -> Qhat
    aux::Dict{String,Float64} # other stuff, among them free energy
end

# constructors
SymmetricPerceptronTSQuenchedGibbsfixedR() = SymmetricPerceptronTSQuenchedGibbsfixedR(1.0,1.0,1.0,1.0,ones(1)*1e-2 ) # dummy constructor for testing
SymmetricPerceptronTSQuenchedGibbsfixedR(β,α,k,R) = SymmetricPerceptronTSQuenchedGibbsfixedR(β,α,k,R,[R])
SymmetricPerceptronTSQuenchedGibbsfixedR(β,α,k,R,O) = SymmetricPerceptronTSQuenchedGibbsfixedR((β=β,α=α,k=k,R=R),NamedArray(O,[:Q]),NamedArray(zeros(1),[:Q]), Dict{String, Float64}())

create_model(Model::SymmetricPerceptronTSQuenchedGibbsfixedR, params) = SymmetricPerceptronTSQuenchedGibbsfixedR(params, Model.O, Model.Oconj, Dict{String, Float64}())

function aux_variables!(Model::SymmetricPerceptronTSQuenchedGibbsfixedR)
    Model.aux["QmR2"] = Model.O[:Q] - Model.params[:R]^2
    Model.aux["1mQ"] = 1 - Model.O[:Q]
    # Model.aux["Rconj"] = 0.1
end

term_A_plus(Model::SymmetricPerceptronTSQuenchedGibbsfixedR, t) = (Model.params[:k] * sqrt(Model.O[:Q]) + Model.params[:R] * t) / sqrt(2 * Model.aux["QmR2"])
term_A_minus(Model::SymmetricPerceptronTSQuenchedGibbsfixedR, t) = (Model.params[:k] * sqrt(Model.O[:Q]) - Model.params[:R] * t) / sqrt(2 * Model.aux["QmR2"])

term_B_plus(Model::SymmetricPerceptronTSQuenchedGibbsfixedR, t) = (Model.params[:k] + sqrt(Model.O[:Q]) * t) / sqrt(2 * Model.aux["1mQ"])
term_B_minus(Model::SymmetricPerceptronTSQuenchedGibbsfixedR, t) = (Model.params[:k] - sqrt(Model.O[:Q]) * t) / sqrt(2 * Model.aux["1mQ"])

num_A_plus(Model::SymmetricPerceptronTSQuenchedGibbsfixedR, t) = (Model.params[:k] * sqrt(Model.O[:Q]) + Model.params[:R] * t)
num_A_minus(Model::SymmetricPerceptronTSQuenchedGibbsfixedR, t) = (Model.params[:k] * sqrt(Model.O[:Q]) - Model.params[:R] * t)

num_A_plus_R(Model::SymmetricPerceptronTSQuenchedGibbsfixedR, t) = Model.params[:R] * (Model.params[:k] * sqrt(Model.O[:Q]) + Model.params[:R] * t)
num_A_minus_R(Model::SymmetricPerceptronTSQuenchedGibbsfixedR, t) = Model.params[:R] * (Model.params[:k] * sqrt(Model.O[:Q]) - Model.params[:R] * t)

num_B_plus(Model::SymmetricPerceptronTSQuenchedGibbsfixedR, t) = (Model.params[:k] + sqrt(Model.O[:Q]) * t)
num_B_minus(Model::SymmetricPerceptronTSQuenchedGibbsfixedR, t) = (Model.params[:k] - sqrt(Model.O[:Q]) * t)

# Conjugate order params SP
function exprOconjall!(Model::SymmetricPerceptronTSQuenchedGibbsfixedR, X::TI) where TI <: IntegrationMethod
    @extract Model params O Oconj aux
    aux_variables!(Model)

    diff_A_q(t) = (1 / sqrt(pi)) * (
    exp(-((term_A_minus(Model, t))^2)) * (Model.params[:k] / (2 * sqrt(2) * sqrt(Model.O[:Q]) *                           sqrt(Model.aux["QmR2"])) - num_A_minus(Model, t) / (2 * sqrt(2) * (Model.aux["QmR2"])^(3/2))) + 
    exp(-((term_A_plus(Model, t))^2)) * (Model.params[:k] / (2 * sqrt(2) * sqrt(Model.O[:Q]) *                             sqrt(Model.aux["QmR2"])) - num_A_plus(Model, t) / (2 * sqrt(2) * (Model.aux["QmR2"])^(3/2))))

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
    
    aux["Rconj"]  = params[:α] * integrand_R
    Oconj[:Q] = - 2 * params[:α] * integrand_q
end

# Order params SP
function exprR(Model::SymmetricPerceptronTSQuenchedGibbsfixedR, X::TI) where TI <: IntegrationMethod
    @extract Model : params Oconj O aux 
    
    IR(z) = tanh(sqrt(Oconj[:Q]) * z + aux["Rconj"])
    return integrate(X, z -> IR(z))
end

function exprQ(Model::SymmetricPerceptronTSQuenchedGibbsfixedR, X::TI) where TI <: IntegrationMethod
    @extract Model : params Oconj aux
    
    IQ(z) = tanh(sqrt(Oconj[:Q]) * z + aux["Rconj"]) ^ 2
    return integrate(X, z -> IQ(z))
end

function lhs!(Model::SymmetricPerceptronTSQuenchedGibbsfixedR, X::TI, Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model :O Oconj
    exprOconjall!(Model, X)
    
    # Onew[:R] = exprR(Model, X)
    eps = 1e-6
    Onew[:Q] = exprQ(Model,X) # clamp(real(exprQ(Model,X)), Onew[:R]^2 + eps, 1 - eps)  
end

########################################################## FREE ENERGY #########################################################


function free_energy(Model::SymmetricPerceptronTSQuenchedGibbsfixedR,X::TI) where TI <: IntegrationMethod
    @extract Model params O Oconj aux
    exprOconjall!(Model,X)

    diff_A_q(t) = (1 / sqrt(pi)) * (
    exp(-((term_A_minus(Model, t))^2)) * (Model.params[:k] / (2 * sqrt(2) * sqrt(Model.O[:Q]) *                           sqrt(Model.aux["QmR2"])) - num_A_minus(Model, t) / (2 * sqrt(2) * (Model.aux["QmR2"])^(3/2))) + 
    exp(-((term_A_plus(Model, t))^2)) * (Model.params[:k] / (2 * sqrt(2) * sqrt(Model.O[:Q]) *                             sqrt(Model.aux["QmR2"])) - num_A_plus(Model, t) / (2 * sqrt(2) * (Model.aux["QmR2"])^(3/2))))

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
    
    # Oconj[:R]  = params[:α] * integrand_R
    aux["Rconj"]  = params[:α] * integrand_R
    Oconj[:Q] = - 2 * params[:α] * integrand_q
    
    
    

    constraints = 0.5 * Oconj[:Q] * (O[:Q] - 1) - Model.params[:R] * aux["Rconj"]
    
    I_GS(t) = log(2 * cosh(sqrt(Oconj[:Q])*t + aux["Rconj"] ))
    
    int1(t) = erfc_A(t) * log(part1_q(t))
    int2(t) = erf_A(t) * log(part2_q(t))
    
    GE = integrate(X,int1) + integrate(X,int2)
    
    G = constraints + integrate(X,I_GS)  + params[:α] * GE
    
    return - G # this quantity is beta*f
end