
# Tricritical (triple) point P / SG / sR of the coupled-students corrupted model, from the linearisation
# of the RS equations around the paramagnet P (m = q = 0, t = t_P != 0). See Tricriticalpoint_corrupted_full.tex.
#
# "Order parameters" updated by lhs!:  βc, αc.
# Auxiliary (Oconj):                   T  = t_P   (pair correlation in P, solved exactly inside lhs! at the current βc)
#                                      Th = t-hat (= αc βc^2 r^4 t/(Δ Δy) at the fixed point)
#
# The two conditions are
#   P -> sR :  λ_m = α βT β r^2 χ / (ΔT Δy) = 1                        (gives αc)
#   P -> SG :  α β^2 μ_q = 1                                          (gives βc, once αc is eliminated)
# with χ = 1 + (y-1) t and
#   general = false : μ_q = χ^2 Cq / y^2  (ansatz (B11) of the paper, the one of Corrupted.jl)
#   general = true  : μ_q = largest eigenvalue of X K D^-2 K (self/cross RS ansatz with private noise)
# Eliminating α with the first condition, the second reads  Λ(β) = β ΔT Δy μ_q/(βT r^2 χ) = 1, with t = t_P(β).
# lhs! solves it exactly (solve_betac!: scan from high T + bisection) and returns (βc, αc), so the outer
# iteration converges at the first step whatever the damping. Diagnostics (λ_m, λ_q, λ_t, κ) in compute_aux!.

struct CoupledHopfieldTS_TCpoint <: FPModel
    params::NamedTuple{(:r,:βT,:γ,:y,:ν),Tuple{FT,FT,FT,Int, Bool}}
    O::NamedVec #vector of orderparameters -> βc, αc
    Oconj::NamedVec #auxiliary -> T (t in P), Th (t-hat in P)
    InnerCW::HomogeneousSpin1Half
    aux::Dict{String,FT} # other stuff
    general::Bool # false: ansatz (B11) (consistent with Corrupted.jl); true: self/cross RS ansatz
end

# constructors
CoupledHopfieldTS_TCpoint() = CoupledHopfieldTS_TCpoint(1.0,0.8,1.0, 2,false,[0.5, 0.1]) # dummy constructor for testing
CoupledHopfieldTS_TCpoint(r,βT,γ,y,ν; general::Bool = false) = CoupledHopfieldTS_TCpoint(r,βT,γ,y,ν, [0.5, 0.1]; general = general)
CoupledHopfieldTS_TCpoint(r,βT,γ,y,ν,O; general::Bool = false) = CoupledHopfieldTS_TCpoint((r=r,βT=βT,γ=γ,y=y,ν=ν),NamedArray(O,[:βc,:αc]),NamedArray([0.0;0.0],[:T,:Th]), HomogeneousSpin1Half(y), Dict{String, FT}(), general)
########################################################## SELF-CONSISTENT EQUATIONS #########################################################

create_model(Model::CoupledHopfieldTS_TCpoint, params) = CoupledHopfieldTS_TCpoint(params,Model.O, Model.Oconj, HomogeneousSpin1Half(Int(params[:y])),Dict{String,Float64}(), Model.general)

set_integrationmethod(Model::CoupledHopfieldTS_TCpoint,X::TI) where TI <: IntegrationMethod = TracedSpin1HalfIntegral(Model,X)

function sanity_checks!(Model::CoupledHopfieldTS_TCpoint)
    Model.params[:βT] >= 1 && throw("You should use a βT < 1 !!!!")
    Model.params[:y] < 2 && throw("Use y >= 2 (for y = 1: αc = (ΔT/(βT r^2))^2, Tc = 1 + ΔT/(βT r^2))")
end

########################################################## SELF-CONSISTENT EQUATIONS #########################################################
DeltaT(Model::CoupledHopfieldTS_TCpoint) =  1 - Model.params[:βT]
Delta(Model::CoupledHopfieldTS_TCpoint, t) =  1 - Model.O[:βc] * ( 1 -  Model.params[:r]^2 * t )
Deltay(Model::CoupledHopfieldTS_TCpoint, t) =  1 - Model.O[:βc] * (1 + Model.params[:r]^2 *(Model.params[:y] - 1) * t )
Delta(Model::CoupledHopfieldTS_TCpoint) = Delta(Model, Model.Oconj[:T])
Deltay(Model::CoupledHopfieldTS_TCpoint) = Deltay(Model, Model.Oconj[:T])

is_free_energy_definite(Model::CoupledHopfieldTS_TCpoint) = (Delta(Model)>0) && (Deltay(Model)>0)

function Deltas!(Model::CoupledHopfieldTS_TCpoint)
    @extract Model : params O Oconj aux
    y, r, t = params[:y], params[:r], Oconj[:T]

    aux["ΔT"] = DeltaT(Model)
    aux["Δ"] = clamp(Delta(Model), 1e-30, Inf)
    aux["Δy"] = clamp(Deltay(Model), 1e-30, Inf)
    aux["chi"] = 1 + (y-1) * t        # χ_∥ in P
    aux["chi_perp"] = 1 - t           # χ_⊥ in P

    aux["Cq"] = ( (y-1) * (1-r^2)^2 / aux["Δ"]^2 ) + ( (1 + r^2 * (y-1) ) ^2 / aux["Δy"]^2 )
end

# zero-field pair correlation τ_y(J) = Ω(ξ1 ξ2) of the y-spin Curie-Weiss block, and its slope
# κ = dτ_y/dJ = (Ω(S^4) - Ω(S^2)^2) / (2y(y-1)) = 1 + 2(y-2)t + (y-2)(y-3)/2 Ω(ξ1ξ2ξ3ξ4) - y(y-1)/2 t^2.
# Max-log-weight shift, so it is safe for any sign of J.
function tau_kappa(C::HomogeneousSpin1Half, J)
    y = C.N
    logw = [0.5 * J * (M^2 - y) + C.logdos[iM] for (iM, M) in enumerate(C.Magnetization)]
    c = maximum(logw)
    Z = 0.0; S2 = 0.0; S4 = 0.0
    for (iM, M) in enumerate(C.Magnetization)
        b = exp(logw[iM] - c)
        Z += b; S2 += b * M^2; S4 += b * M^4
    end
    S2 /= Z; S4 /= Z
    return (S2 - y) / (y * (y - 1)), (S4 - S2^2) / (2 * y * (y - 1))
end

# J_eff in P on the P-sR line, where α is eliminated: t-hat = α β^2 r^4 t/(Δ Δy) = β r^4 t ΔT/(βT r^2 χ Δ)
function Jeff_sR(Model::CoupledHopfieldTS_TCpoint, t)
    @extract Model : params O
    β, r, y = O[:βc], params[:r], params[:y]
    beta_nu_prefactor = params[:ν] ? β : 1.0
    that = β * r^2 * t * DeltaT(Model) / ( params[:βT] * (1 + (y-1)*t) * Delta(Model, t) )
    return beta_nu_prefactor * params[:γ] / y + that, that
end

# t in P at the current βc: smallest root of t = τ_y(J_eff(t)), the branch connected to high temperature.
# The map is increasing for Δy > 0, so iterating from t = 0 converges monotonically to it.
function expr_OconjT!(Model::CoupledHopfieldTS_TCpoint; tol = 1e-14, maxiter = 100000)
    @extract Model : Oconj
    t = 0.0
    for it in 1:maxiter
        tnew = tau_kappa(Model.InnerCW, Jeff_sR(Model, t)[1])[1]
        if abs(tnew - t) < tol
            t = tnew
            break
        end
        t = tnew
    end
    Oconj[:T] = t
    Oconj[:Th] = Jeff_sR(Model, t)[2]
    Deltas!(Model)
end

# μ_q: SG eigenvalue per unit α β^2 in P (used for the diagnostics).
#   (B11): (1,1) entry  χ^2 Cq / y^2.
#   general: largest eigenvalue of X K D^-2 K, X = diag(χ_∥^2, χ_⊥^2), D^-2 = diag(Δy^-2, Δ^-2),
#            K = (1/y) [1+(y-1)r^2  (y-1)(1-r^2); 1-r^2  y-1+r^2],  det K = r^2.
function expr_muq(Model::CoupledHopfieldTS_TCpoint)
    @extract Model : params aux
    y, r = params[:y], params[:r]
    χp, χo, Δ, Δy = aux["chi"], aux["chi_perp"], aux["Δ"], aux["Δy"]
    a11 = χp^2 * aux["Cq"] / y^2
    Model.general || return a11
    a22 = χo^2 * ( (y-1) * (1-r^2)^2 / Δy^2 + (y-1+r^2)^2 / Δ^2 ) / y^2
    detT = (χp * χo * r^2)^2 / (Δ * Δy)^2
    tr = a11 + a22
    return 0.5 * ( tr + sqrt(max(tr^2 - 4 * detT, 0.0)) )
end

# SG eigenvalue on the P-sR line: Λ(β) = α_sR(β) β^2 μ_q with α_sR = ΔT Δy/(βT r^2 χ β), t = t_P(β).
# Λ -> 0 for β -> 0 and Λ -> ∞ for Δy -> 0+. Sets O[:βc] = β and the auxiliaries at that β.
function Lambda_sR!(Model::CoupledHopfieldTS_TCpoint, β)
    @extract Model : params O aux
    O[:βc] = β
    expr_OconjT!(Model)
    (Delta(Model) > 0 && Deltay(Model) > 0) || return Inf
    return β * aux["ΔT"] * aux["Δy"] * expr_muq(Model) / ( params[:βT] * params[:r]^2 * aux["chi"] )
end

# P-SG condition Λ(β) = 1 with α eliminated by the P-sR one.
# (B11) it reads  β = βT r^2 y^2/(ΔT Δy χ Cq), but iterating that form is unstable (slope ≈ -3 at y = 2,
# since Δy is small at the tricritical point). Here: first crossing of Λ(β) = 1 coming from high
# temperature (scan on a grid, then bisection). Returns NaN if there is no crossing with Δy > 0.
function solve_betac!(Model::CoupledHopfieldTS_TCpoint; ngrid = 400, tol = 1e-15)
    lo, hi = 0.0, NaN
    for β in range(1e-3, 1 - 1e-9, length = ngrid)
        if Lambda_sR!(Model, β) >= 1
            hi = β
            break
        end
        lo = β
    end
    isnan(hi) && return NaN
    while hi - lo > tol
        mid = 0.5 * (lo + hi)
        Lambda_sR!(Model, mid) >= 1 ? (hi = mid) : (lo = mid)
    end
    β = 0.5 * (lo + hi)
    Lambda_sR!(Model, β)        # leave t, Δ, Δy, χ evaluated at βc
    return β
end

# P-sR condition
function expr_alphac(Model::CoupledHopfieldTS_TCpoint)
    @extract Model : params O aux
    return aux["ΔT"] * aux["Δy"] / ( params[:βT] * params[:r]^2 *  aux["chi"] * O[:βc] )
end

# Onew = exact (βc, αc) for the current parameters, so the outer iteration converges at once.
function lhs!(Model::CoupledHopfieldTS_TCpoint,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    β_old = Model.O[:βc]
    βc = solve_betac!(Model)
    isnan(βc) && @warn "no tricritical point with Δy > 0 for params = $(Model.params)"
    Onew[:βc] = βc
    Onew[:αc] = expr_alphac(Model)     # at O[:βc] = βc
    Model.O[:βc] = β_old
end

# Diagnostics at the current (βc, αc): at the tricritical point λ_m = λ_q = 1 and λ_t < 1.
function compute_aux!(Model::CoupledHopfieldTS_TCpoint,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj aux
    expr_OconjT!(Model)
    β, α, r, y, t = O[:βc], O[:αc], params[:r], params[:y], Oconj[:T]
    J = Jeff_sR(Model, t)[1]
    κ = tau_kappa(Model.InnerCW, J)[2]
    aux["Tc"] = 1 / β
    aux["t"] = t
    aux["Jeff"] = J
    aux["kappa"] = κ
    aux["lambda_m"] = α * params[:βT] * β * r^2 * aux["chi"] / (aux["ΔT"] * aux["Δy"])    # P -> sR
    aux["lambda_q"] = α * β^2 * expr_muq(Model)                                          # P -> SG
    aux["lambda_t"] = κ * α * β^2 * r^4 * ( (y-1) / aux["Δy"]^2 + 1 / aux["Δ"]^2 ) / y     # stability of t_P
end

########################################################## FREE ENERGY #########################################################

function free_energy(Model::CoupledHopfieldTS_TCpoint,X::TI) where TI <: IntegrationMethod
    return -1.
end
