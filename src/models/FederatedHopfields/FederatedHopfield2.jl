struct FederatedHopfield2 <: FPModel
    params::NamedTuple{(:β,:α1,:α2,:γ,:βT,:ν),Tuple{FT,FT,FT,FT,FT,Int}}
    O::NamedVec #vector of orderparameters -> M
    Oconj::NamedVec #vector of orderparameters -> M
    aux::Dict{String,FT} # other stuff, among them free energy
end

#
FederatedHopfield2() = FederatedHopfield2(.8,0.1,.1,.0,0.8,0,ones(4)*0.3) # dummy constructor for testing
FederatedHopfield2(β,α1,α2,γ,βT,ν) = FederatedHopfield2(β,α1,α2,γ,βT,ν,ones(4)) 
FederatedHopfield2(β,α1,α2,γ,βT,ν,O) = FederatedHopfield2((β=β,α1=α1,α2=α2,γ=γ,βT=βT,ν=ν),NamedArray(O,[:M1,:M2,:Q1,:Q2]),NamedArray(zeros(4),[:M1,:M2,:Q1,:Q2]), Dict{String, FT}())

create_model(Model::FederatedHopfield2, params) = FederatedHopfield2(params, Model.O,Model.Oconj,Dict{String, FT}())

########################################################## SELF-CONSISTENT EQUATIONS #########################################################

function aux_variables!(Model::FederatedHopfield2)
    @extract Model : O Oconj params aux 
    aux["ΔT"] = 1 - params[:βT]
    aux["Δ1"] = 1 - params[:β] * (1-O[:Q1])
    aux["Δ2"] = 1 - params[:β] * (1-O[:Q2])
end


function exprOconjall!(Model::FederatedHopfield2)
    @extract Model : O Oconj params aux
    aux_variables!(Model)

    Oconj[:M1] = ( params[:α1] * params[:βT] * params[:β] * O[:M1] )  / (aux["ΔT"] * aux["Δ1"] )
    Oconj[:M2] = ( params[:α2] * params[:βT] * params[:β] * O[:M2] )  / (aux["ΔT"] * aux["Δ2"] )

    
    Oconj[:Q1] = ( params[:α1] * params[:βT] * (params[:β]^2) * (O[:M1]^2) )  / (aux["ΔT"] * aux["Δ1"]^2 ) + params[:α1] * params[:β]^2 * O[:Q1]  /  aux["Δ1"]^2

    Oconj[:Q2] = ( params[:α2] * params[:βT] * (params[:β]^2) * (O[:M2]^2) )  / (aux["ΔT"] * aux["Δ2"]^2 ) + params[:α2] * params[:β]^2 * O[:Q2]  /  aux["Δ2"]^2

    
end


# order params SP
function internal_mag1(Model::FederatedHopfield2,z1,z2)
    @extract Model : params Oconj aux
    return tanh( Oconj[:M1] + sqrt(Oconj[:Q1])*z1  + atanh(tanh(0.5 * params[:γ] * params[:β]^params[:ν]) * tanh(Oconj[:M2] + sqrt(Oconj[:Q2])*z2) ) )
end


function internal_mag2(Model::FederatedHopfield2,z1,z2)
    @extract Model : params Oconj aux
    return tanh( Oconj[:M2] + sqrt(Oconj[:Q2])*z2  + atanh(tanh(0.5 * params[:γ] * params[:β]^params[:ν]) * tanh(Oconj[:M1] + sqrt(Oconj[:Q1])*z1) ) )
end

function exprM1(Model::FederatedHopfield2,X::TI) where TI <: IntegrationMethod
    IM1(z1,z2) = internal_mag1(Model,z1,z2)
    return integrate_double(X,IM1)
end

function exprM2(Model::FederatedHopfield2,X::TI) where TI <: IntegrationMethod
    @extract Model : params Oconj aux
    IM2(z1,z2) = internal_mag2(Model,z1,z2)
    return integrate_double(X,IM2)
end



function exprQ1(Model::FederatedHopfield2,X::TI) where TI <: IntegrationMethod
    IQ1(z1,z2) = internal_mag1(Model,z1,z2)^2
    return integrate_double(X,IQ1)
end

function exprQ2(Model::FederatedHopfield2,X::TI) where TI <: IntegrationMethod
    @extract Model : params Oconj aux
    IQ2(z1,z2) = internal_mag2(Model,z1,z2)^2
    return integrate_double(X,IQ2)
end


function lhs!(Model::FederatedHopfield2,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O
    exprOconjall!(Model)
    Onew[:M1] = exprM1(Model,X)
    Onew[:M2] = exprM2(Model,X)
    Onew[:Q1] = exprQ1(Model,X)
    Onew[:Q2] = exprQ2(Model,X)
end

########################################################## FREE ENERGY #########################################################
function free_energy(Model::FederatedHopfield2,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj aux 
    exprOconjall!(Model)

    piece_overlap = (0.5/params[:β])*Oconj[:Q1] * (1-O[:Q1]) + (0.5/params[:β])*Oconj[:Q2] * (1-O[:Q2])
    piece_magnetizations = O[:M1] * Oconj[:M1] / params[:β] + O[:M2] * Oconj[:M2] / params[:β]
    log_det = 0.5 * (params[:α1] / params[:β]) * (log(aux["ΔT"]) + log(aux["Δ1"])) + 0.5 * (params[:α2] / params[:β]) * (log(aux["ΔT"]) + log(aux["Δ2"])) 

    other_1 = 0.5 * params[:α1] * O[:Q1] / aux["Δ1"] + 0.5 * params[:α1] * (params[:βT] * O[:M1]^2) / (aux["ΔT"] * aux["Δ1"])

    other_2 = 0.5 * params[:α2] * O[:Q2] / aux["Δ2"] + 0.5 * params[:α2] * (params[:βT] * O[:M2]^2) / (aux["ΔT"] * aux["Δ2"])

    entropy = 

    Ientropy(z1,z2) = log(sum( exp( ξ1*( Oconj[:M1] + sqrt(Oconj[:Q1])*z1 ) + ξ2*( Oconj[:M2] + sqrt(Oconj[:Q2])*z2 )  + 0.5 * params[:γ] * params[:β]^params[:ν]* ξ1 * ξ2) for ξ1 in [-1;1] for ξ2 in [-1;1]))

    entropy = integrate_double(X,Ientropy)
    return piece_overlap + piece_magnetizations + log_det - other_1 - other_2 - entropy/params[:β]

end


_isSignal(Model::FederatedHopfield2, X::TI, Oold::NamedVec,tol::FT) where TI <: IntegrationMethod = (Model.O[:M1] > 1e-1) & (Model.O[:M2] > 1e-1)

########################################################## RS (dAT) STABILITY #########################################################
# Replicon (dAT) stability of the RS solution of the SDD scenario (App. C of arXiv:2606.13045): two students with separate
# datasets of loads α1, α2. Derivation, checks and notation: dAT_SDN/dAT_SDD.pdf (script: dAT_SDN/verify_dAT_SDD.py).
#
# The energetic term -Σ_u (α_u/2) ln det Ξ_u depends only on the self-overlaps q^u_ab, so the replicon modes are
# δq^u_ab = ε^ab v_u and the stability operator acting on (v1, v2) is the 2x2 matrix
#     T = G_S G_E,   G_E = diag(g1, g2),  g_u = α_u β² / Δ_u²,  Δ_u = 1 - β(1 - Q_u),
#                    G_S = [A1 B; B A2],  A_u = E_z[χ_uu²],  B = E_z[χ_12²],
# with χ the connected correlation matrix of the two-spin block P(ξ1,ξ2) ∝ exp(h1 ξ1 + h2 ξ2 + J ξ1 ξ2),
# h_u = M̂_u + sqrt(Q̂_u) z_u (z1, z2 independent N(0,1)), J = γ β^ν / 2. RS is stable iff
#     Λ+ = [g1 A1 + g2 A2 + sqrt((g1 A1 - g2 A2)² + 4 g1 g2 B²)] / 2 < 1,
# and the dAT line is Λ+ = 1, i.e. (Δ1² - α1 β² A1)(Δ2² - α2 β² A2) = α1 α2 β⁴ B² with Δ_u² > α_u β² A_u.
# At γ = 0, B = 0 and the condition splits into the AGS conditions of the two students.

lse_sdd(a, b) = max(a, b) + log1p(exp(-abs(a - b)))          # log(e^a + e^b)

# (χ11, χ22, χ12) of the two-spin block, closed form evaluated in log domain (no overflow at large fields):
# χ11 = 2[cosh 2J + cosh 2h2]/Z², χ22 = 2[cosh 2J + cosh 2h1]/Z², χ12 = 2 sinh(2J)/Z², Z = e^J cosh(h1+h2) + e^-J cosh(h1-h2)
function twospin_chis(h1, h2, J)
    l2 = log(2) - 2 * lse_sdd(J + logcosh(h1 + h2), -J + logcosh(h1 - h2))       # log(2/Z²)
    χ11 = exp(l2 + lse_sdd(logcosh(2J), logcosh(2h2)))
    χ22 = exp(l2 + lse_sdd(logcosh(2J), logcosh(2h1)))
    χ12 = sign(J) * exp(l2 - log(2) + 2abs(J) + log1p(-exp(-4abs(J))))           # 2|sinh 2J| = e^{2|J|}(1 - e^{-4|J|})
    return χ11, χ22, χ12
end

"""
    dat_sdd2!(aux, X, β, α1, α2, J, Oconj, Δ1, Δ2)

2x2 replicon operator of two SDD students, shared by FederatedHopfield2 and FederatedHopfield2SameAlpha (α1 = α2).
Oconj = (M1, M2, Q1, Q2) conjugates at the current point, J = γ β^ν / 2. Stored in aux:
- "A1", "A2", "B": entropic averages E[χ11²], E[χ22²], E[χ12²]
- "Lambda_1" = g1 A1, "Lambda_2" = g2 A2: diagonal entries of T (each student's AGS eigenvalue if B were zero)
- "Lambda_dAT" = Λ+ (largest eigenvalue of T), "Lambda_dAT_minus" = Λ-
  (on a student-symmetric solution: Λ± = g (A ± B), eigenvectors (1, ±1))
- "dAT" = (Δ1² - α1β²A1)(Δ2² - α2β²A2) - α1α2β⁴B²: division-free margin, changes sign on the dAT line
- "dAT_v2_over_v1": unstable direction δq²/δq¹ = (Λ+ - g1A1)/(g2 B)
- "RS_valid" (Δ1, Δ2 > 0), "is_RS_stable" (RS_valid and Λ+ < 1)
NaN-safe: at a non-finite point every Λ is NaN and is_RS_stable = false.
"""
function dat_sdd2!(aux, X::TI, β, α1, α2, J, Oconj, Δ1, Δ2) where TI <: IntegrationMethod
    h1(z) = Oconj[:M1] + sqrt(Oconj[:Q1]) * z
    h2(z) = Oconj[:M2] + sqrt(Oconj[:Q2]) * z
    A1 = integrate_double(X, (z1, z2) -> twospin_chis(h1(z1), h2(z2), J)[1]^2)
    A2 = integrate_double(X, (z1, z2) -> twospin_chis(h1(z1), h2(z2), J)[2]^2)
    B  = integrate_double(X, (z1, z2) -> twospin_chis(h1(z1), h2(z2), J)[3]^2)

    g1, g2 = α1 * β^2 / Δ1^2, α2 * β^2 / Δ2^2
    L1, L2 = g1 * A1, g2 * A2
    disc = sqrt((L1 - L2)^2 + 4 * g1 * g2 * B^2)

    aux["A1"], aux["A2"], aux["B"] = A1, A2, B
    aux["Lambda_1"], aux["Lambda_2"] = L1, L2
    aux["Lambda_dAT"] = 0.5 * (L1 + L2 + disc)
    aux["Lambda_dAT_minus"] = 0.5 * (L1 + L2 - disc)
    aux["dAT"] = (Δ1^2 - α1 * β^2 * A1) * (Δ2^2 - α2 * β^2 * A2) - α1 * α2 * β^4 * B^2
    aux["dAT_v2_over_v1"] = g2 * B > 0 ? (aux["Lambda_dAT"] - L1) / (g2 * B) : 0.0
    aux["RS_valid"] = (Δ1 > 0) && (Δ2 > 0)
    aux["is_RS_stable"] = (aux["RS_valid"] == 1.0) && aux["Lambda_dAT"] < 1.0
    return nothing
end

# Replicon (dAT) stability of the current RS point (O = M1, M2, Q1, Q2); aux keys as in dat_sdd2!
function eigs_dat!(Model::FederatedHopfield2, X::TI) where TI <: IntegrationMethod
    @extract Model : params Oconj aux
    exprOconjall!(Model)                                    # Oconj and Δ1, Δ2 at the current O
    β = params[:β]
    J = 0.5 * params[:γ] * β^params[:ν]
    dat_sdd2!(aux, X, β, params[:α1], params[:α2], J, Oconj, aux["Δ1"], aux["Δ2"])
end

# called by solve_SPequations! at convergence (after free_energy!)
compute_aux!(Model::FederatedHopfield2, X::TI) where TI <: IntegrationMethod = eigs_dat!(Model, X)

# keepgoing function for critical_line.jl: scanning a parameter stops when RS becomes unstable (dAT line)
_is_RS_stable(Model::FederatedHopfield2, X::TI, Oold::NamedVec, tol::FT) where TI <: IntegrationMethod = Model.aux["Lambda_dAT"] < 1.0

# P -> SG line (replicon = longitudinal instability of the paramagnet, A_u = 1, B = tanh²J). With x = (T-1)²:
# (x - α1)(x - α2) = α1 α2 tanh⁴J, largest root. ν = 0: closed form; ν ≠ 0: J depends on T, solved by bisection in
# T ∈ [1 + sqrt(max α), 1 + sqrt(α1 + α2)], which always brackets the root since tanh⁴J ∈ [0, 1].
function T_PSG(α1, α2, γ, ν = 0)
    xroot(C) = 0.5 * (α1 + α2 + sqrt((α1 - α2)^2 + 4 * α1 * α2 * C^4))
    ν == 0 && return 1 + sqrt(xroot(tanh(0.5 * γ)))
    lo, hi = 1 + sqrt(max(α1, α2)), 1 + sqrt(α1 + α2)
    lo == hi && return lo
    return find_zero(T -> (1 + sqrt(xroot(tanh(0.5 * γ * T^(-ν))))) - T, (lo, hi), Roots.Bisection())
end
T_PSG(Model::FederatedHopfield2) = T_PSG(Model.params[:α1], Model.params[:α2], Model.params[:γ], Model.params[:ν])
