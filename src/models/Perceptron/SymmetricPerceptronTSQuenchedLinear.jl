struct SymmetricPerceptronTSQuenchedLinear <: FPModel
    params::NamedTuple{(:β,:α,:k),Tuple{FT,FT,FT}}
    O::NamedVec #vector of orderparameters -> M, Q
    Oconj::NamedVec #vector of conjugate orderparameters -> Qhat
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors
SymmetricPerceptronTSQuenchedLinear() = SymmetricPerceptronTSQuenchedLinear(1.0,1.0,1.0,ones(2)*1e-2 ) # dummy constructor for testing
SymmetricPerceptronTSQuenchedLinear(β,α,k) = SymmetricPerceptronTSQuenchedLinear(β,α,k,zeros(2))
SymmetricPerceptronTSQuenchedLinear(β,α,k,O) = SymmetricPerceptronTSQuenchedLinear((β=β,α=α,k=k),NamedArray(O,[:R;:Q]),NamedArray(zeros(2),[:R;:Q]), Dict{String, FT}())

create_model(Model::SymmetricPerceptronTSQuenchedLinear, params) = SymmetricPerceptronTSQuenchedLinear(params, Model.O, Model.Oconj, Dict{String, FT}())

function aux_variables!(Model::SymmetricPerceptronTSQuenchedLinear)
    eps = 1e-36
    Model.aux["QmR2"] = clamp(Model.O[:Q] - Model.O[:R]^2, eps, Inf)
    Model.aux["1mQ"] = clamp(1 - Model.O[:Q], eps, Inf)
end

term_A_plus(Model::SymmetricPerceptronTSQuenchedLinear, t) = (Model.params[:k] * sqrt(Model.O[:Q]) + (Model.O[:R]) * t) / sqrt(2 * Model.aux["QmR2"])
term_A_minus(Model::SymmetricPerceptronTSQuenchedLinear, t) = (Model.params[:k] * sqrt(Model.O[:Q]) - (Model.O[:R]) * t) / sqrt(2 * Model.aux["QmR2"])

term_B_plus(Model::SymmetricPerceptronTSQuenchedLinear, t) = (Model.params[:k] + sqrt(Model.O[:Q]) * t) / sqrt(2 * Model.aux["1mQ"])
term_B_minus(Model::SymmetricPerceptronTSQuenchedLinear, t) = (Model.params[:k] - sqrt(Model.O[:Q]) * t) / sqrt(2 * Model.aux["1mQ"])

num_A_plus(Model::SymmetricPerceptronTSQuenchedLinear, t) = (Model.params[:k] * sqrt(Model.O[:Q]) + (Model.O[:R]) * t)
num_A_minus(Model::SymmetricPerceptronTSQuenchedLinear, t) = (Model.params[:k] * sqrt(Model.O[:Q]) - (Model.O[:R]) * t)

num_A_plus_R(Model::SymmetricPerceptronTSQuenchedLinear, t) = Model.O[:R] * (Model.params[:k] * sqrt(Model.O[:Q]) + (Model.O[:R]) * t)
num_A_minus_R(Model::SymmetricPerceptronTSQuenchedLinear, t) = Model.O[:R] * (Model.params[:k] * sqrt(Model.O[:Q]) - (Model.O[:R]) * t)

num_B_plus(Model::SymmetricPerceptronTSQuenchedLinear, t) = (Model.params[:k] + sqrt(Model.O[:Q]) * t)
num_B_minus(Model::SymmetricPerceptronTSQuenchedLinear, t) = (Model.params[:k] - sqrt(Model.O[:Q]) * t)

num_B_k_plus(Model::SymmetricPerceptronTSQuenchedLinear,k, t) = (k + sqrt(Model.O[:Q]) * t + Model.params[:β]*(1-Model.O[:Q])) / sqrt(2 * Model.aux["1mQ"])
num_B_k_minus(Model::SymmetricPerceptronTSQuenchedLinear,k, t) = (k + sqrt(Model.O[:Q]) * t - Model.params[:β]*(1-Model.O[:Q])) / sqrt(2 * Model.aux["1mQ"])

# Conjugate order params SP
function exprOconjall!(Model::SymmetricPerceptronTSQuenchedLinear, X::TI) where TI <: IntegrationMethod
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

    diff_B_k_minus(k,t) = (sqrt(Model.O[:Q]) * Model.aux["1mQ"]*(Model.params[:β])+t+k*sqrt(Model.O[:Q]) ) / (2 * sqrt(Model.O[:Q]) * (Model.aux["1mQ"])^(3/2))
    diff_B_k_plus(k,t) = (sqrt(Model.O[:Q]) * Model.aux["1mQ"]*(-Model.params[:β])+t+k*sqrt(Model.O[:Q]) ) / (2 * sqrt(Model.O[:Q]) * (Model.aux["1mQ"])^(3/2))
    
    P1_(t) = 0.5*exp(-Model.params[:β]*sqrt(Model.O[:Q])*t)*(erf(num_B_k_minus(Model,Model.params[:k],t)) - erf(num_B_k_minus(Model,0.0,t))) + 
        0.5*exp(Model.params[:β]*sqrt(Model.O[:Q])*t)*(erf(num_B_k_plus(Model,0.0,t)) - erf(num_B_k_plus(Model,-Model.params[:k],t)))

    P1(t) = 0.5*exp(-Model.params[:β]*sqrt(Model.O[:Q])*t)*(erf(num_B_k_minus(Model,0.0,t), num_B_k_minus(Model,Model.params[:k],t)))  +
          0.5*exp(Model.params[:β]*sqrt(Model.O[:Q])*t)*(erf(num_B_k_plus(Model,-Model.params[:k],t), num_B_k_plus(Model,0.0,t)) )

    P2_(t) = cosh(Model.params[:β]*sqrt(Model.O[:Q])*t) - 0.5*exp(Model.params[:β]*sqrt(Model.O[:Q])*t)*erf(num_B_k_plus(Model,Model.params[:k],t)) + 
        0.5*exp(-Model.params[:β]*sqrt(Model.O[:Q])*t)*erf(num_B_k_minus(Model,-Model.params[:k],t))

    P2(t) =  0.5*exp(Model.params[:β]*sqrt(Model.O[:Q])*t)*erfc(num_B_k_plus(Model,Model.params[:k],t)) +
             0.5*exp(-Model.params[:β]*sqrt(Model.O[:Q])*t)*erfc(-num_B_k_minus(Model,-Model.params[:k],t))

    #diff_P1_Q(t) = - 0.5*Model.params[:β]*t/sqrt(2*pi*Model.O[:Q])*exp(-Model.params[:β]*sqrt(Model.O[:Q])*t)*(exp(-num_B_k_minus(Model,Model.params[:k],t)^2)*diff_B_k_minus(Model.params[:k],t) - 
    #                exp(-num_B_k_minus(Model,0.0,t)^2)*diff_B_k_minus(0.0,t) )  +
    #               0.5*Model.params[:β]*t/sqrt(2*pi*Model.O[:Q])*exp(Model.params[:β]*sqrt(Model.O[:Q])*t)*(exp(-num_B_k_plus(Model,Model.params[:k],t)^2)*diff_B_k_plus(Model.params[:k],t) - 
    #                exp(-num_B_k_plus(Model,0.0,t)^2)*diff_B_k_plus(0.0,t) ) 

    diff_P1_Q(t) =  -(Model.params[:β]*t)/(4*sqrt(Model.O[:Q]))*exp(-Model.params[:β]*sqrt(Model.O[:Q])*t)*(erf(num_B_k_minus(Model,0.0,t), num_B_k_minus(Model,Model.params[:k],t)) ) +
                    (Model.params[:β]*t)/(4*sqrt(Model.O[:Q]))*exp(Model.params[:β]*sqrt(Model.O[:Q])*t)*(erf(num_B_k_plus(Model,-Model.params[:k],t), num_B_k_plus(Model,0.0,t)) ) +
                    1.0/sqrt(2*pi)*exp(-Model.params[:β]*sqrt(Model.O[:Q])*t)*(exp(-num_B_k_minus(Model,Model.params[:k],t)^2)*diff_B_k_minus(Model.params[:k],t) - 
                    exp(-num_B_k_minus(Model,0.0,t)^2)*diff_B_k_minus(0.0,t) )  +
                    1.0/sqrt(2*pi)*exp(Model.params[:β]*sqrt(Model.O[:Q])*t)*(exp(-num_B_k_plus(Model,0.0,t)^2)*diff_B_k_plus(0.0,t) - 
                    exp(-num_B_k_plus(Model,-Model.params[:k],t)^2)*diff_B_k_plus(-Model.params[:k],t) )                    

    #diff_P2_Q(t) = 0.5*Model.params[:β]*t/sqrt(Model.O[:Q])*sinh(Model.params[:β]*sqrt(Model.O[:Q])*t) +
    #    0.5*Model.params[:β]*t/sqrt(2*pi*Model.O[:Q])*exp(-Model.params[:β]*sqrt(Model.O[:Q])*t)*exp(-num_B_k_plus(Model,Model.params[:k],t)^2)*diff_B_k_plus(Model.params[:k],t) +
    #    0.5*Model.params[:β]*t/sqrt(2*pi*Model.O[:Q])*exp(Model.params[:β]*sqrt(Model.O[:Q])*t)*exp(-num_B_k_minus(Model,-Model.params[:k],t)^2)*diff_B_k_minus(-Model.params[:k],t)
    
    diff_P2_Q_(t) = 0.5*Model.params[:β]*t/sqrt(Model.O[:Q])*sinh(Model.params[:β]*sqrt(Model.O[:Q])*t) -
         Model.params[:β]*t/(4*sqrt(Model.O[:Q]))*exp(Model.params[:β]*sqrt(Model.O[:Q])*t)*erf(num_B_k_plus(Model,Model.params[:k],t)) -
         Model.params[:β]*t/(4*sqrt(Model.O[:Q]))*exp(-Model.params[:β]*sqrt(Model.O[:Q])*t)*erf(num_B_k_minus(Model,-Model.params[:k],t)) -
         1.0/sqrt(2*pi)*exp(Model.params[:β]*sqrt(Model.O[:Q])*t)*exp(-num_B_k_plus(Model,Model.params[:k],t)^2)*diff_B_k_plus(Model.params[:k],t) +
         1.0/sqrt(2*pi)*exp(-Model.params[:β]*sqrt(Model.O[:Q])*t)*exp(-num_B_k_minus(Model,-Model.params[:k],t)^2)*diff_B_k_minus(-Model.params[:k],t)
    
    diff_P2_Q(t) = 0.25*Model.params[:β]*t/sqrt(Model.O[:Q])*(erfc(num_B_k_plus(Model,Model.params[:k],t))*exp(Model.params[:β]*sqrt(Model.O[:Q])*t) -
                                                              erfc(-num_B_k_minus(Model,-Model.params[:k],t))*exp(-Model.params[:β]*sqrt(Model.O[:Q])*t) ) -
         1.0/sqrt(2*pi)*exp(Model.params[:β]*sqrt(Model.O[:Q])*t)*exp(-num_B_k_plus(Model,Model.params[:k],t)^2)*diff_B_k_plus(Model.params[:k],t) +
         1.0/sqrt(2*pi)*exp(-Model.params[:β]*sqrt(Model.O[:Q])*t)*exp(-num_B_k_minus(Model,-Model.params[:k],t)^2)*diff_B_k_minus(-Model.params[:k],t)

    erf_A(t) = 0.5 * (erf(term_A_plus(Model, t)) + erf(term_A_minus(Model, t)))
    
    erfc_A(t) = 0.5 * (1 - erf(term_A_plus(Model, t))) + 0.5 * (1 - erf(term_A_minus(Model, t)))
    
    erf_B(t) = 0.5 * (erf(term_B_plus(Model, t)) + erf(term_B_minus(Model, t)))
    
    erfc_B(t) = 0.5 * (1 - erf(term_B_plus(Model, t))) + 0.5 * (1 - erf(term_B_minus(Model, t)))
    
    
    A1(t) = exp(-Model.params[:β]*Model.params[:k]+0.5*(1-Model.O[:Q])*Model.params[:β]^2) * P1(t)
    A1_(t) = exp(-Model.params[:β]*Model.params[:k]+0.5*(1-Model.O[:Q])*Model.params[:β]^2) * P1_(t)
    A2(t) = exp(Model.params[:β]*Model.params[:k]+0.5*(1-Model.O[:Q])*Model.params[:β]^2) * P2(t)
    A2_(t) = exp(Model.params[:β]*Model.params[:k]+0.5*(1-Model.O[:Q])*Model.params[:β]^2) * P2_(t)

    B1(t) = exp(-Model.params[:β]) * erf_B(t)
    B2(t) = exp(-Model.params[:β]) * erfc_B(t)
    #println(integrate(X,A1)," ",integrate(X,A1_)," ",integrate(X,A2)," ",integrate(X,A2_))
    #println(P1(1)," ",P1(2)," ",P1(5)," ",P1(10))
    #println(P1_(1)," ",P1_(2)," ",P1_(5)," ",P1_(10))
    #println(P2(1)," ",P2(2)," ",P2(5)," ",P2(10))
    #println(P2_(1)," ",P2_(2)," ",P2_(5)," ",P2_(10))

    #println(integrate(X,B1)," ",integrate(X,B2))
    #println(Model.O)

    DP2 = integrate(X,diff_P2_Q)
    DP2_ = integrate(X,diff_P2_Q_)
    #println(DP2," ",DP2_)

    part1_q(t) = erfc_B(t) + exp(-Model.params[:β]*Model.params[:k]+0.5*(1-Model.O[:Q])*Model.params[:β]^2) * P1(t)
    # part1_q(t) = erfc_B(t) + exp(-Model.params[:β]) * erf_B(t)
    # part2_q(t) = exp(-Model.params[:β]) * erfc_B(t) + erf_B(t)
    part2_q(t) = exp(Model.params[:β]*Model.params[:k]+0.5*(1-Model.O[:Q])*Model.params[:β]^2) * P2(t) + erf_B(t)

    part1_q_V0(t) = erfc_B(t) + exp(-Model.params[:β]) * erf_B(t)
    part2_q_V0(t) = exp(-Model.params[:β]) * erfc_B(t) + erf_B(t)
    f1_q_V0(t) = - diff_A_q(t) * log(part1_q_V0(t))
    f2_q_V0(t) = erfc_A(t) * ((- diff_B_q(t) + exp(-Model.params[:β]) * diff_B_q(t)) / part1_q_V0(t))
    f3_q_V0(t) = diff_A_q(t) * log(part2_q_V0(t))
    f4_q_V0(t) = erf_A(t) * ((- exp(-Model.params[:β]) * diff_B_q(t) + diff_B_q(t)) / part2_q_V0(t))

    
    f1_q(t) = - diff_A_q(t) * log(part1_q(t))
    
    # f2_q(t) = erfc_A(t) * ((- diff_B_q(t) + exp(-Model.params[:β]) * diff_B_q(t)) / part1_q(t))
    f2_q(t) = erfc_A(t) * ((- diff_B_q(t) + exp(-Model.params[:β]*Model.params[:k]+0.5*(1-Model.O[:Q])*Model.params[:β]^2)*( diff_P1_Q(t) - P1(t)*0.5*Model.params[:β]^2)) / part1_q(t))
    
    f3_q(t) = diff_A_q(t) * log(part2_q(t))
    
    f4_q(t) = erf_A(t) * ((exp(Model.params[:β]*Model.params[:k]+0.5*(1-Model.O[:Q])*Model.params[:β]^2)*(diff_P2_Q(t) - P2(t)*0.5*Model.params[:β]^2) + diff_B_q(t)) / part2_q(t))

    integrand_q = integrate(X,f1_q) + integrate(X,f2_q) + integrate(X,f3_q) + integrate(X,f4_q)
    #println(integrand_q)
    #println(integrate(X,f1_q)," ",integrate(X,f2_q)," ",integrate(X,f3_q)," ",integrate(X,f4_q))
    #println(integrate(X,f1_q_V0)," ",integrate(X,f2_q_V0)," ",integrate(X,f3_q_V0)," ",integrate(X,f4_q_V0))
    #println()    

    f1_R(t) = - diff_A_R(t) * log(part1_q(t))
    
    f2_R(t) = diff_A_R(t) * log(part2_q(t))
    
    integrand_R = integrate(X,f1_R) + integrate(X,f2_R)
    
    Oconj[:R]  = params[:α] * integrand_R
    eps = 1e-15
    Oconj[:Q] = clamp(- 2 * params[:α] * integrand_q, eps, Inf)
    # println(Oconj[:Q])
    #  Oconj[:Q] = clamp(- 2 * params[:α] * integrand_q, eps, Inf)
end

# Order params SP
function exprR(Model::SymmetricPerceptronTSQuenchedLinear, X::TI) where TI <: IntegrationMethod
    @extract Model : params Oconj O 
    O[:R] == 0.0 && return 0.0
    
    IR(z) = tanh(sqrt(Oconj[:Q]) * z + Oconj[:R])
    return integrate(X, z -> IR(z))
end

function exprQ(Model::SymmetricPerceptronTSQuenchedLinear, X::TI) where TI <: IntegrationMethod
    @extract Model : params Oconj 
    
    IQ(z) = tanh(sqrt(Oconj[:Q]) * z + Oconj[:R]) ^ 2
    return integrate(X, z -> IQ(z))
end

function lhs!(Model::SymmetricPerceptronTSQuenchedLinear, X::TI, Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model :O Oconj
    exprOconjall!(Model, X)
    
    Onew[:R] = exprR(Model, X)
    Onew[:Q] = exprQ(Model,X) #clamp(real(exprQ(Model,X)), Onew[:R]^2 + eps, 1 - eps)  
end

########################################################## FREE ENERGY #########################################################


function free_energy(Model::SymmetricPerceptronTSQuenchedLinear,X::TI) where TI <: IntegrationMethod
    @extract Model params O Oconj aux
    exprOconjall!(Model,X)




    P1(t) = 0.5*exp(-Model.params[:β]*sqrt(Model.O[:Q])*t)*(erf(num_B_k_minus(Model,0.0,t), num_B_k_minus(Model,Model.params[:k],t)))  +
          0.5*exp(Model.params[:β]*sqrt(Model.O[:Q])*t)*(erf(num_B_k_plus(Model,-Model.params[:k],t), num_B_k_plus(Model,0.0,t)) )

    P2(t) =  0.5*exp(Model.params[:β]*sqrt(Model.O[:Q])*t)*erfc(num_B_k_plus(Model,Model.params[:k],t)) +
             0.5*exp(-Model.params[:β]*sqrt(Model.O[:Q])*t)*erfc(-num_B_k_minus(Model,-Model.params[:k],t))

    erf_A(t) = 0.5 * (erf(term_A_plus(Model, t)) + erf(term_A_minus(Model, t)))
    
    erfc_A(t) = 0.5 * (1 - erf(term_A_plus(Model, t))) + 0.5 * (1 - erf(term_A_minus(Model, t)))
    
    erf_B(t) = 0.5 * (erf(term_B_plus(Model, t)) + erf(term_B_minus(Model, t)))
    
    erfc_B(t) = 0.5 * (1 - erf(term_B_plus(Model, t))) + 0.5 * (1 - erf(term_B_minus(Model, t)))
    
    # part1_q(t) = erfc_B(t) + exp(-Model.params[:β]) * erf_B(t)
    part1_q(t) = erfc_B(t) + exp(-Model.params[:β]*Model.params[:k]+0.5*(1-Model.O[:Q])*Model.params[:β]^2) * P1(t)
    part2_q(t) = exp(Model.params[:β]*Model.params[:k]+0.5*(1-Model.O[:Q])*Model.params[:β]^2) * P2(t) + erf_B(t)
    
    #P1_FE(t) = erfc_A(t)*log(part1_q(t))
    #P2_FE(t) = erf_A(t)*log(part2_q(t))

    #Oconj[:R]  = params[:α] * integrand_R
    #eps = 1e-15
    #Oconj[:Q] = clamp(- 2 * params[:α] * integrand_q, eps, Inf)
    

    constraints = 0.5 * Oconj[:Q] * (O[:Q] - 1) - O[:R] * Oconj[:R]
    #println(Model.O)
    #println(Model.Oconj)
    I_GS(t) = log(2 * cosh(sqrt(Oconj[:Q])*t + Oconj[:R] ))
    
    int1(t) = erfc_A(t) * log(part1_q(t))
    int2(t) = erf_A(t) * log(part2_q(t))
    
    GE = integrate(X,int1) + integrate(X,int2)
    
    G = constraints + integrate(X,I_GS)  + params[:α] * GE
    
    return - G # this quantity is beta*f
end
