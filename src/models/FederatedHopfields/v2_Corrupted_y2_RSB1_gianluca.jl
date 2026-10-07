# =============================================================================
# Model with the Uniform Corruption (UC)
# =============================================================================

# struct v2_Corrupted_2S_1RSB <: FPModel
#     params::NamedTuple{
#         (:β, :α, :βT, :r1, :r2, :γ, :ν),
#         Tuple{Float64, Float64, Float64, Float64, Float64, Float64, Int}
#     }
#     O::NamedVec
#     Oconj::NamedVec
#     aux::Dict{String,Float64}
# end

# # -----------------------------------------------------------------------------
# # Constructors
# # -----------------------------------------------------------------------------

# v2_Corrupted_2S_1RSB() =
#     v2_Corrupted_2S_1RSB(0.1, 0.1, 0.1, 1.0, 1.0, 2.0, 2, ones(6))

# v2_Corrupted_2S_1RSB(β, α, βT, r1, r2, γ, ν) =
#     v2_Corrupted_2S_1RSB(β, α, βT, r1, r2, γ, ν, zeros(6))

# v2_Corrupted_2S_1RSB(β, α, βT, r1, r2, γ, ν, O) =
#     v2_Corrupted_2S_1RSB(
#         (β=β, α=α, βT=βT, r1=r1, r2=r2, γ=γ, ν=ν),
#         NamedArray(O, [:m1, :m2, :q1, :q2, :q, :t]),
#         NamedArray(zeros(6), [:m1, :m2, :q1, :q2, :q, :t]),
#         Dict{String,Float64}()
#     )

# create_model(Model::v2_Corrupted_2S_1RSB, params) =
#     v2_Corrupted_2S_1RSB(params, Model.O, Model.Oconj, Dict{String,Float64}())

# # =============================================================================
# # Self-consistent equations
# # =============================================================================

# DeltaT(Model::v2_Corrupted_2S_1RSB) = 1 - Model.params[:βT]

# function Deltas!(Model::v2_Corrupted_2S_1RSB)
#     @extract Model : O aux params

#     aux["ΔT"] = DeltaT(Model)

#     aux["det_log"] =
#         (1 - params[:β] * (1 - O[:q1])) *
#         (1 - params[:β] * (1 - O[:q2])) -
#         (params[:r1] * params[:r2] * params[:β] * (O[:t] - O[:q]))^2

#     # SAFE boolean check (no Float64 used as Bool)
#     # if aux["det_log"] < 0
#     #     println("det_log < 0")
#     # end

    
# end

# # -----------------------------------------------------------------------------
# # Conjugate order parameters
# # -----------------------------------------------------------------------------

# function exprOconj!(Model::v2_Corrupted_2S_1RSB)
#     @extract Model : O Oconj aux params
#     Deltas!(Model)

#     ΔT = aux["ΔT"]

#     denom =
#         ΔT + params[:β] * ΔT * (
#             -2 + O[:q1] + O[:q2] +
#             params[:β] * (
#                 1 + O[:q1] * (-1 + O[:q2]) - O[:q2] -
#                 params[:r1]^2 * params[:r2]^2 * (O[:q] - O[:t])^2
#             )
#         )

#     Oconj[:m1] =
#         params[:α] * params[:β] * params[:βT] * params[:r1]^2 *
#         (O[:m1] +
#          params[:β] * O[:m1] * (-1 + O[:q2]) +
#          params[:β] * O[:m2] * params[:r2]^2 * (-O[:q] + O[:t])) / denom

#     Oconj[:m2] =
#         params[:α] * params[:β] * params[:βT] * params[:r2]^2 *
#         (O[:m2] +
#          params[:β] * O[:m2] * (-1 + O[:q1]) +
#          params[:β] * O[:m1] * params[:r1]^2 * (-O[:q] + O[:t])) / denom

#     # ---- q1 ----
#     Oconj[:q1] =
#         params[:α] * params[:β] *
#         (-ΔT * (1 + params[:β] * (-1 + 2 * O[:q2])) -
#          params[:β] * params[:βT] * O[:m2]^2 * params[:r2]^2) / denom

#     # ---- q2 ----
#     Oconj[:q2] =
#         -params[:α] * params[:β]^2 *
#         (-ΔT * (1 + params[:β] * (-1 + O[:q1]))^2 * O[:q2]) /
#         (ΔT * (1 + params[:β] * (-2 + O[:q1] + O[:q2]))^2)

#     # ---- q ----
#     Oconj[:q] =
#         params[:α] * params[:β]^2 * params[:r1]^2 * params[:r2]^2 *
#         (ΔT * O[:q]) /
#         (ΔT * (1 + params[:β] * (-2 + O[:q1] + O[:q2]))^2)

#     # ---- t ----
#     Oconj[:t] =
#         params[:α] * params[:β]^2 * params[:r1]^2 * params[:r2]^2 *
#         (ΔT * O[:t]) /
#         (ΔT * (1 + params[:β] * (-2 + O[:q1] + O[:q2]))^2)

    
# end

# =============================================================================
# Manual inner trace
# =============================================================================

# function inner_trace_manual(Model; nbins=500, bound=20.0)
#     # Extract parameters
#     Oconj = Model.Oconj
#     params = Model.params

#     J = Oconj[:t] - Oconj[:q] + params[:γ] / 2
#     det_cov = Oconj[:q1] * Oconj[:q2] - Oconj[:q]^2
#     sqrt_det_cov = sqrt(abs(det_cov))  

#     # Uniform grid and step size
#     Avals = range(-bound, bound; length=nbins)
#     dA = step(Avals)   # spacing

#     xi_vals = (-1, 1)

#     # Initialize accumulators
#     norm = 0.0
#     M1 = 0.0
#     M2 = 0.0
#     Q1 = 0.0
#     Q2 = 0.0
#     Q = 0.0
#     T = 0.0

#     # Loop over 2D grid
#     for A1 in Avals, A2 in Avals
#         # Gaussian weight
#         exponent = - (A1^2 * Oconj[:q2] - 2 * A1 * A2 * Oconj[:q] + A2^2 * Oconj[:q1]) / (2 * det_cov)
#         weight = (exp(exponent) )/(sqrt_det_cov) # include grid cell area no dA

#         # Local fields
#         h1 = Oconj[:m1] + A1
#         h2 = Oconj[:m2] + A2

#         # Inner sum over xi1, xi2 = ±1
#         Z = 0.0
#         m1_inner = 0.0
#         m2_inner = 0.0
#         chi_inner = 0.0

#         for xi1 in xi_vals, xi2 in xi_vals
#             w = exp(h1 * xi1 + h2 * xi2 + J * xi1 * xi2)
#             Z += w
#             m1_inner += xi1 * w
#             m2_inner += xi2 * w
#             chi_inner += xi1 * xi2 * w
#         end

#         # Normalize inner averages
#         m1_inner /= Z
#         m2_inner /= Z
#         chi_inner /= Z

#         # Accumulate weighted sums
#         norm += weight
#         M1 += m1_inner * weight
#         M2 += m2_inner * weight
#         Q1 += (m1_inner^2) * weight
#         Q2 += (m2_inner^2) * weight
#         Q  += (m1_inner * m2_inner) * weight
#         T  += chi_inner * weight
#     end

#     # println("Sanity check: total norm ≈ 2π√detΣ = ", 2*pi*sqrt(det_cov), ", numeric norm = ", norm)

#     # Return normalized expectation values
#     return M1 / norm, M2 / norm, Q1 / norm, Q2 / norm, Q / norm, T / norm
# end

# function inner_trace_diag(Model::v2_Corrupted_2S_1RSB, X::LegendreQuadrature)
#     O = Model.Oconj
#     p = Model.params

#     # ------------------------------------------------------------------
#     # Covariance eigen-decomposition
#     # ------------------------------------------------------------------
#     q1 = O[:q1]
#     q2 = O[:q2]
#     q  = O[:q]

#     Δ = sqrt((q1 - q2)^2 + 4*q^2)
#     λplus  = 0.5 * (q1 + q2 + Δ)
#     λminus = 0.5 * (q1 + q2 - Δ)

#     θ = 0.5 * atan(2*q / (q1 - q2))
#     c = cos(θ)
#     s = sin(θ)

#     lp = sqrt(λplus)
#     lm = sqrt(λminus)

#     # Interaction
#     J = O[:t] - O[:q] + p[:γ]/2

#     # Quadrature fields
#     zs = X.points
#     ws = X.weights
#     SW = X.sumweights

#     # Accumulators
#     M1 = 0.0; M2 = 0.0
#     Q1 = 0.0; Q2 = 0.0
#     Q  = 0.0; T  = 0.0
#     norm = 0.0

#     # ------------------------------------------------------------------
#     # Gauss–Legendre double integral
#     # ------------------------------------------------------------------
#     for (z1, w1) in zip(zs, ws)
#         for (z2, w2) in zip(zs, ws)

#             # combined weight for 2D integral
#             weight = (w1 * w2) / SW
#             norm += weight  # accumulate normalization

#             # transform to A1, A2
#             A1 = c*lp*z1 - s*lm*z2
#             A2 = s*lp*z1 + c*lm*z2

#             h1 = O[:m1] + A1
#             h2 = O[:m2] + A2

#             # local spin trace
#             Z = 0.0
#             m1 = 0.0
#             m2 = 0.0
#             chi = 0.0

#             for xi1 in (-1,1), xi2 in (-1,1)
#                 w = exp(h1*xi1 + h2*xi2 + J*xi1*xi2)
#                 Z += w
#                 m1 += xi1*w
#                 m2 += xi2*w
#                 chi += xi1*xi2*w
#             end

#             m1 /= Z
#             m2 /= Z
#             chi /= Z

#             # accumulate order parameters
#             M1 += weight*m1
#             M2 += weight*m2
#             Q1 += weight*m1^2
#             Q2 += weight*m2^2
#             Q  += weight*m1*m2
#             T  += weight*chi
#         end
#     end

#     return M1, M2, Q1, Q2, Q, T
# end
##############################################

# @inline function cholesky2x2(q1, q2, q)
#     println(q1)
#     q1 > 0 || error("q1 ≤ 0 in cholesky2x2")
#     l11 = sqrt(q1)
#     l21 = q / l11
#     l22sq = q2 - l21^2
#     l22sq > 0 || error("Covariance matrix not positive definite")
#     println("END cholesky2x2")
#     return l11, l21, sqrt(l22sq)
# end

# function inner_trace_GH(Model::v2_Corrupted_2S_1RSB; nGH=20)
#     Oconj  = Model.Oconj
#     params = Model.params

#     q1 = Oconj[:q1]
#     q2 = Oconj[:q2]
#     q  = Oconj[:q]

#     # ---- diagonalize covariance ----
#     l11, l21, l22 = cholesky2x2(q1, q2, q)

#     # ---- GH nodes and weights ----
#     x, w = gausshermite(nGH)

#     xi_vals = (-1, 1)
#     J = Oconj[:t] - Oconj[:q] + params[:γ] / 2

#     norm = 0.0
#     M1 = 0.0
#     M2 = 0.0
#     Q1 = 0.0
#     Q2 = 0.0
#     Q  = 0.0
#     T  = 0.0

#     for i in 1:nGH, j in 1:nGH
#         # standard normal variables
#         z1 = sqrt(2) * x[i]
#         z2 = sqrt(2) * x[j]
#         weight = w[i] * w[j] / Base.MathConstants.pi   # correct Dz normalization

#         # correlated Gaussian fields
#         A1 = l11 * z1
#         A2 = l21 * z1 + l22 * z2

#         h1 = Oconj[:m1] + A1
#         h2 = Oconj[:m2] + A2

#         Z = 0.0
#         m1loc = 0.0
#         m2loc = 0.0
#         chiloc = 0.0

#         for ξ1 in xi_vals, ξ2 in xi_vals
#             wξ = exp(h1*ξ1 + h2*ξ2 + J*ξ1*ξ2)
#             Z += wξ
#             m1loc += ξ1 * wξ
#             m2loc += ξ2 * wξ
#             chiloc += ξ1 * ξ2 * wξ
#         end

#         m1loc /= Z
#         m2loc /= Z
#         chiloc /= Z

#         norm += weight
#         M1 += m1loc * weight
#         M2 += m2loc * weight
#         Q1 += m1loc^2 * weight
#         Q2 += m2loc^2 * weight
#         Q  += m1loc*m2loc * weight
#         T  += chiloc * weight
#     end

#     return M1 / norm,M2 / norm,Q1 / norm,  Q2 / norm, Q  / norm,  T  / norm
# end


# # =============================================================================
# # Left-hand side update
# # =============================================================================

# function lhs!(Model::v2_Corrupted_2S_1RSB, X::TI, Onew::NamedVec) where TI <: IntegrationMethod
#     exprOconj!(Model)
#     @assert iszero(Onew)

#     vals = inner_trace_GH(Model )
#     Onew[:m1], Onew[:m2], Onew[:q1], Onew[:q2], Onew[:q], Onew[:t] = real.(vals)

# end


#########2dotentativo
struct v2_Corrupted_2S_1RSB <: FPModel
    params::NamedTuple{(:β,:α,:βT,:r1,:r2,:γ,:ν),Tuple{Float64,Float64,Float64,Float64,Float64,Float64, Int}}
    O::NamedVec #vector of orderparameters -> m1, m2, q1,q2,q, t
    Oconj::NamedVec #vector of conjugate orderparameters -> Qha
    aux::Dict{String,Float64} # other stuff, among them free energy
end

# constructors
v2_Corrupted_2S_1RSB() = v2_Corrupted_2S_1RSB(.1, .1,.1,1.0,1.0,2, 2, ones(6)) # dummy constructor for testing
v2_Corrupted_2S_1RSB(β,α,βT,r1,r2,γ,ν) = v2_Corrupted_2S_1RSB(β,α,βT,r1,r2,γ,ν,zeros(6))

# In the next line I change HomogeneousSpin1Half(y)---->deleted
v2_Corrupted_2S_1RSB(β,α,βT,r1,r2,γ,ν,O) = v2_Corrupted_2S_1RSB((β=β,α=α,βT=βT,r1=r1,r2=r2,γ=γ,ν=ν),NamedArray(O,[:m1;:m2;:q1;:q2;:q;:t]),NamedArray(zeros(6),[:m1;:m2;:q1;:q2;:q;:t]),Dict{String,Float64}())

create_model(Model::v2_Corrupted_2S_1RSB, params) = v2_Corrupted_2S_1RSB(params,Model.O, Model.Oconj, Dict{String,Float64}())


########################################################## SELF-CONSISTENT EQUATIONS #########################################################
###### USEFUL COMBINATION
### DENOMINATORS of SPEQs
DeltaT(Model::v2_Corrupted_2S_1RSB) =  1 - Model.params[:βT]

function Deltas!(Model::v2_Corrupted_2S_1RSB)
    @extract Model : O Oconj aux params
    @assert params[:ν] == 0
    Model.aux["ΔT"] = DeltaT(Model)

    Model.aux["det_log"] = (1-params[:β]*(1-O[:q1])) * (1-params[:β]*(1-O[:q2])) - (params[:r1]*params[:r2]*params[:β] * (O[:t] - O[:q]))^2

    # Model.aux["det_log"] < 0 && println("det_log <0")
end



#### Conjugate order params 
function exprOconj!(Model::v2_Corrupted_2S_1RSB)
    @extract Model : O Oconj aux params
    Deltas!(Model)
    
    Oconj[:m1] =  (params[:α] * params[:β] * params[:βT] * params[:r1]^2 * 
                    (O[:m1] + params[:β] * O[:m1] * (-1 + O[:q2]) + 
                    params[:β] * O[:m2] * params[:r2]^2 * (-O[:q] + O[:t]))) / 
                    (aux["ΔT"] + 
                    params[:β] * aux["ΔT"] * 
                    (-2 + O[:q1] + O[:q2] + 
                    params[:β] * (1 + O[:q1] * (-1 + O[:q2]) - O[:q2] - 
                                    params[:r1]^2 * params[:r2]^2 * (O[:q] - O[:t])^2)))
    
    
    Oconj[:m2] = (
    params[:α] * params[:β] * params[:βT] * params[:r2]^2 *
    (
        O[:m2] +
        params[:β] * O[:m2] * (-1 + O[:q1]) +
        params[:β] * O[:m1] * params[:r1]^2 * (-O[:q] + O[:t])
    )
) / (
    aux["ΔT"] +
    params[:β] * aux["ΔT"] * (
        -2 + O[:q1] + O[:q2] +
        params[:β] * (
            1 + O[:q1] * (-1 + O[:q2]) - O[:q2] -
            params[:r1]^2 * params[:r2]^2 * (O[:q] - O[:t])^2
        )
    )
)
    
    Oconj[:q1] =  params[:α] * (
    (
        params[:β] * (
            -aux["ΔT"] * (1 + params[:β] * (-1 + 2 * O[:q2])) -
            params[:β] * params[:βT] * O[:m2]^2 * params[:r2]^2
        )
    ) / (
        aux["ΔT"] + params[:β] * aux["ΔT"] * (
            -2 + O[:q1] + O[:q2] +
            params[:β] * (
                1 + O[:q1] * (-1 + O[:q2]) - O[:q2] -
                params[:r1]^2 * params[:r2]^2 * (O[:q] - O[:t])^2
            )
        )
    ) +
    (
        params[:β] * (1 + params[:β] * (-1 + O[:q2]))
    ) / (
        (1 + params[:β] * (-1 + O[:q1])) *
        (1 + params[:β] * (-1 + O[:q2])) -
        (params[:β] * O[:q] * params[:r1] * params[:r2] -
         params[:β] * params[:r1] * params[:r2] * O[:t])^2
    ) -
    (
        params[:β]^2 * aux["ΔT"] * (1 + params[:β] * (-1 + O[:q2])) *
        (
            -aux["ΔT"] * (
                O[:q1] + O[:q2] + params[:β] * O[:q1] * (-1 + 2 * O[:q2]) -
                params[:β] * (
                    O[:q2] + 2 * O[:q] * params[:r1]^2 * params[:r2]^2 * (O[:q] - O[:t])
                )
            ) -
            params[:βT] * (
                O[:m1]^2 * (1 + params[:β] * (-1 + O[:q2])) * params[:r1]^2 +
                O[:m2]^2 * (1 + params[:β] * (-1 + O[:q1])) * params[:r2]^2 +
                2 * params[:β] * O[:m1] * O[:m2] * params[:r1]^2 * params[:r2]^2 * (-O[:q] + O[:t])
            )
        )
    ) / (
        aux["ΔT"] + params[:β] * aux["ΔT"] * (
            -2 + O[:q1] + O[:q2] +
            params[:β] * (
                1 + O[:q1] * (-1 + O[:q2]) - O[:q2] -
                params[:r1]^2 * params[:r2]^2 * (O[:q] - O[:t])^2
            )
        )
    )^2
)
    
    # (-params[:α])*final_derivative(:q1, O[:m1], O[:m2], O[:q1], O[:q2], O[:t], O[:q], params[:β], params[:βT], aux["ΔT"], params[:r1], params[:r2])

    Oconj[:q2] =  -(
                    params[:α] * params[:β]^2 * (
                        -aux["ΔT"] * (1 + params[:β] * (-1 + O[:q1]))^2 * O[:q2] +
                        params[:β] * aux["ΔT"] * params[:r1]^2 * params[:r2]^2 * (O[:q] - O[:t]) *
                        (
                            O[:q] * (2 + params[:β] * (-2 + O[:q1])) +
                            params[:β] * O[:q1] * O[:t]
                        ) -
                        params[:βT] * params[:r2]^2 *
                        (
                            O[:m2] + params[:β] * O[:m2] * (-1 + O[:q1]) +
                            params[:β] * O[:m1] * params[:r1]^2 * (-O[:q] + O[:t])
                        )^2
                    )
                ) / (
                    aux["ΔT"] * (
                        1 + params[:β] * (
                            -2 + O[:q1] + O[:q2] +
                            params[:β] * (
                                1 + O[:q1] * (-1 + O[:q2]) - O[:q2] -
                                params[:r1]^2 * params[:r2]^2 * (O[:q] - O[:t])^2
                            )
                        )
                    )^2
                )
    Oconj[:q] =  (
                    params[:α] * params[:β]^2 * params[:r1]^2 * params[:r2]^2 * (
                        aux["ΔT"] * O[:q] -
                        params[:βT] * (
                            O[:m1] * (-1 + params[:β] - params[:β] * O[:q2]) +
                            params[:β] * O[:m2] * params[:r2]^2 * (O[:q] - O[:t])
                        ) * (
                            O[:m2] + params[:β] * O[:m2] * (-1 + O[:q1]) +
                            params[:β] * O[:m1] * params[:r1]^2 * (-O[:q] + O[:t])
                        ) +
                        params[:β] * aux["ΔT"] * (
                            params[:β] * O[:q]^3 * params[:r1]^2 * params[:r2]^2 +
                            (O[:q1] + O[:q2] - params[:β] * O[:q2] + params[:β] * O[:q1] * (-1 + 2 * O[:q2])) * O[:t] -
                            2 * params[:β] * O[:q]^2 * params[:r1]^2 * params[:r2]^2 * O[:t] +
                            O[:q] * (
                                -2 + params[:β] - params[:β] * O[:q1] * O[:q2] +
                                params[:β] * params[:r1]^2 * params[:r2]^2 * O[:t]^2
                            )
                        )
                    )
                ) / (
                    aux["ΔT"] * (
                        1 + params[:β] * (
                            -2 + O[:q1] + O[:q2] +
                            params[:β] * (
                                1 + O[:q1] * (-1 + O[:q2]) - O[:q2] -
                                params[:r1]^2 * params[:r2]^2 * (O[:q] - O[:t])^2
                            )
                        )
                    )^2
                )
    Oconj[:t] =  (
                    params[:α] * params[:β]^2 * params[:r1]^2 * params[:r2]^2 * (
                        -params[:βT] * (
                            O[:m1] * (-1 + params[:β] - params[:β] * O[:q2]) +
                            params[:β] * O[:m2] * params[:r2]^2 * (O[:q] - O[:t])
                        ) * (
                            O[:m2] + params[:β] * O[:m2] * (-1 + O[:q1]) +
                            params[:β] * O[:m1] * params[:r1]^2 * (-O[:q] + O[:t])
                        ) +
                        aux["ΔT"] * (
                            -params[:β] * O[:q] * (O[:q1] + O[:q2]) +
                            O[:t] +
                            2 * params[:β] * (-1 + O[:q1] + O[:q2]) * O[:t] +
                            params[:β]^2 * (
                                2 * O[:q]^3 * params[:r1]^2 * params[:r2]^2 +
                                O[:t] -
                                5 * O[:q]^2 * params[:r1]^2 * params[:r2]^2 * O[:t] -
                                O[:t] * (
                                    O[:q1] * (2 - 3 * O[:q2]) +
                                    2 * O[:q2] +
                                    params[:r1]^2 * params[:r2]^2 * O[:t]^2
                                ) +
                                O[:q] * (
                                    O[:q1] + O[:q2] - 2 * O[:q1] * O[:q2] +
                                    4 * params[:r1]^2 * params[:r2]^2 * O[:t]^2
                                )
                            )
                        )
                    )
                ) / (
                    aux["ΔT"] * (
                        1 + params[:β] * (
                            -2 + O[:q1] + O[:q2] +
                            params[:β] * (
                                1 + O[:q1] * (-1 + O[:q2]) - O[:q2] -
                                params[:r1]^2 * params[:r2]^2 * (O[:q] - O[:t])^2
                            )
                        )
                    )^2
                )
end

 # Use same Gauss-Hermite structure but rescale for variance 1/2
function rescaled_for_half_variance(X::LegendreQuadrature)
    points = X.points ./ sqrt(2)
    weights = X.weights .* sqrt(2)
    return LegendreQuadrature(X.n, X.extr, points, weights, sum(weights))
end

# Integral with the 4 variables z1,z2, w, bar w

function inner_trace(Model::v2_Corrupted_2S_1RSB, X::LegendreQuadrature; n_uv=20, bound_uv=5.0)
    @extract Model : params O Oconj aux

    Xu_raw = LegendreQuadrature(n_uv; bound=bound_uv)
    Xv_raw = LegendreQuadrature(n_uv; bound=bound_uv)

    Xu = rescaled_for_half_variance(Xu_raw)
    Xv = rescaled_for_half_variance(Xv_raw)

    # Prepare accumulators (make sure to initialize before loops!)
    M1 = 0.0 + 0.0im
    M2 = 0.0 + 0.0im
    Q1  = 0.0 + 0.0im
    Q2  = 0.0 + 0.0im
    Q  = 0.0 + 0.0im
    T  = 0.0 + 0.0im
    

    J = Oconj[:t] - Oconj[:q] + params[:γ]/2

    # Iterate over all quadrature points (nested 4 loops)
    for (iz1, z1) in enumerate(X.points)
        w1 = X.weights[iz1]
        for (iz2, z2) in enumerate(X.points)
            w2 = X.weights[iz2]

            for (iu, u) in enumerate(Xu.points)
                w3 = Xu.weights[iu]
                for (iv, v) in enumerate(Xv.points)
                    w4 = Xv.weights[iv]

                    # Define complex omega variables
                    ω = u + im*v
                    ω_bar = u - im*v

                    # Compute h1, h2 using given params and quadrature points
                    h1 = Oconj[:m1] + z1 * sqrt(Complex(Oconj[:q1])) + ω * sqrt(Complex(Oconj[:q]))
                    h2 = Oconj[:m2] + z2 * sqrt(Complex(Oconj[:q2])) + ω_bar * sqrt(Complex(Oconj[:q]))
                    
                    # Inner average over discrete xi1, xi2 = ±1
                    xi_vals = [-1, 1]
                    Z = 0.0 + 0.0im
                    m1_inner = 0.0 + 0.0im
                    m2_inner = 0.0 + 0.0im
                    chi_inner = 0.0 + 0.0im

                    for xi1 in xi_vals, xi2 in xi_vals
                        weight_ = exp(h1 * xi1 + h2 * xi2 + J * xi1 * xi2)
                        Z += weight_
                        m1_inner += xi1 * weight_
                        m2_inner += xi2 * weight_
                        chi_inner += xi1 * xi2 * weight_
                    end

                    # Normalize inner averages
                    m1_inner /= Z
                    m2_inner /= Z
                    chi_inner /= Z

                    # Now integrate these quantities weighted by product of weights
                    W = w1 * w2 * w3 * w4

                    M1 += m1_inner * W
                    M2 += m2_inner * W
                    Q1 += m1_inner^2 *W
                    Q2 += m2_inner^2 *W
                    Q += m1_inner*m2_inner *W
                    T  += chi_inner * W

                end
            end
        end
    end

    normalization = X.sumweights^2 * Xu.sumweights * Xv.sumweights

    return M1 / normalization, M2 / normalization, Q1 / normalization, Q2 / normalization, Q / normalization, T / normalization
end

# integral with A1, A2
function inner_trace_manual(Model; nbins=500, bound=10.0)
    Oconj = Model.Oconj
    params = Model.params

    J = Oconj[:t] - Oconj[:q] + params[:γ] / 2
    det_cov = Oconj[:q1] * Oconj[:q2] - Oconj[:q]^2
    if det_cov <= 0
        error("Covariance not positive definite: det_cov=$det_cov")
    end
    sqrt_det_cov = sqrt(abs(det_cov))

    Avals = range(-bound, bound; length=nbins)
    dA = step(Avals)
    xi_vals = (-1, 1)

    norm = M1 = M2 = Q1 = Q2 = Q = T = 0.0

    for A1 in Avals, A2 in Avals
        exponent = - (A1^2 * Oconj[:q2] - 2 * A1 * A2 * Oconj[:q] + A2^2 * Oconj[:q1]) / (2 * det_cov)
        weight = exp(exponent) * (1 / 2*sqrt_det_cov) 

        h1 = Oconj[:m1] + A1
        h2 = Oconj[:m2] + A2

        Z = m1_inner = m2_inner = chi_inner = 0.0
        for xi1 in xi_vals, xi2 in xi_vals
            w = exp(h1 * xi1 + h2 * xi2 + J * xi1 * xi2)
            Z += w
            m1_inner += xi1 * w
            m2_inner += xi2 * w
            chi_inner += xi1 * xi2 * w
        end
        m1_inner /= Z
        m2_inner /= Z
        chi_inner /= Z

        norm += weight
        M1 += m1_inner * weight
        M2 += m2_inner * weight
        Q1 += m1_inner^2 * weight
        Q2 += m2_inner^2 * weight
        Q  += m1_inner * m2_inner * weight
        T  += chi_inner * weight
    end

    return M1 / norm, M2 / norm, Q1 / norm, Q2 / norm, Q / norm, T / norm, norm
end


# function inner_trace_manual_diag(Model; nbins=500, bound=20.0)
#     Oconj = Model.Oconj
#     params = Model.params

#     # Eigen-decomposition of covariance
#     det_cov = Oconj[:q1]*Oconj[:q2] - Oconj[:q]^2
#     λ_plus  = 0.5 * (Oconj[:q1] + Oconj[:q2] + sqrt((Oconj[:q1]-Oconj[:q2])^2 + 4*Oconj[:q]^2))
#     λ_minus = 0.5 * (Oconj[:q1] + Oconj[:q2] - sqrt((Oconj[:q1]-Oconj[:q2])^2 + 4*Oconj[:q]^2))
#     θ = 0.5 * atan(2*Oconj[:q] / (Oconj[:q1] - Oconj[:q2]))

#     # Quadrature points for standard normals
#     z_vals = range(-bound, bound; length=nbins)
#     dz = step(z_vals)

#     xi_vals = (-1, 1)
#     M1 = 0.0; M2 = 0.0; Q1 = 0.0; Q2 = 0.0; Q = 0.0; T = 0.0
#     norm = 0.0

#     J = Oconj[:t] - Oconj[:q] + params[:γ]/2

#     for z1 in z_vals
#         for z2 in z_vals
#             # Diagonalized linear transformation
#             A1 = cos(θ)*sqrt(λ_plus)*z1 - sin(θ)*sqrt(λ_minus)*z2
#             A2 = sin(θ)*sqrt(λ_plus)*z1 + cos(θ)*sqrt(λ_minus)*z2

#             h1 = Oconj[:m1] + A1
#             h2 = Oconj[:m2] + A2

#             # Inner sum over discrete spins
#             Z = 0.0; m1_inner = 0.0; m2_inner = 0.0; chi_inner = 0.0
#             for xi1 in xi_vals, xi2 in xi_vals
#                 w = exp(h1*xi1 + h2*xi2 + J*xi1*xi2)
#                 Z += w
#                 m1_inner += xi1*w
#                 m2_inner += xi2*w
#                 chi_inner += xi1*xi2*w
#             end
#             m1_inner /= Z
#             m2_inner /= Z
#             chi_inner /= Z

#             # Standard Gaussian weight
#             weight = exp(-(z1^2 + z2^2)/2) / (2π)

#             # Riemann sum element
#             W = weight * dz * dz

#             M1 += m1_inner*W
#             M2 += m2_inner*W
#             Q1 += m1_inner^2*W
#             Q2 += m2_inner^2*W
#             Q  += m1_inner*m2_inner*W
#             T  += chi_inner*W
#             norm += W
#         end
#     end

#     return M1/norm, M2/norm, Q1/norm, Q2/norm, Q/norm, T/norm, norm
# end


# function inner_trace_legendre_diag(Model::v2_Corrupted_2S_1RSB,X::TI) where TI <: IntegrationMethod
#     Oconj = Model.Oconj
#     params = Model.params

#     # Eigen-decomposition of covariance
#     λ_plus  = 0.5 * (Oconj[:q1] + Oconj[:q2] + sqrt((Oconj[:q1]-Oconj[:q2])^2 + 4*Oconj[:q]^2))
#     λ_minus = 0.5 * (Oconj[:q1] + Oconj[:q2] - sqrt((Oconj[:q1]-Oconj[:q2])^2 + 4*Oconj[:q]^2))
#     θ = 0.5 * atan(2*Oconj[:q] / (Oconj[:q1] - Oconj[:q2]))

#     xi_vals = (-1, 1)
#     M1 = 0.0; M2 = 0.0; Q1 = 0.0; Q2 = 0.0; Q = 0.0; T = 0.0

#     # Coupling term
#     J = Oconj[:t] - Oconj[:q] + params[:γ]/2

#     # Double integration using Legendre quadrature
#     for (iz1, z1) in enumerate(X.points)
#         w1 = X.weights[iz1]
#         for (iz2, z2) in enumerate(X.points)
#             w2 = X.weights[iz2]

#             # Transform to correlated A1, A2 using rotation
#             A1 = cos(θ)*sqrt(λ_plus)*z1 - sin(θ)*sqrt(λ_minus)*z2
#             A2 = sin(θ)*sqrt(λ_plus)*z1 + cos(θ)*sqrt(λ_minus)*z2

#             h1 = Oconj[:m1] + A1
#             h2 = Oconj[:m2] + A2

#             # Inner sum over discrete spins
#             Z = 0.0
#             m1_inner = 0.0
#             m2_inner = 0.0
#             chi_inner = 0.0

#             for xi1 in xi_vals, xi2 in xi_vals
#                 w = exp(h1*xi1 + h2*xi2 + xi1*xi2*J)
#                 Z += w
#                 m1_inner += xi1*w
#                 m2_inner += xi2*w
#                 chi_inner += xi1*xi2*w
#             end

#             m1_inner /= Z
#             m2_inner /= Z
#             chi_inner /= Z

#             # Accumulate weighted
#             W = w1 * w2
#             M1 += m1_inner * W
#             M2 += m2_inner * W
#             Q1 += m1_inner^2 * W
#             Q2 += m2_inner^2 * W
#             Q  += m1_inner*m2_inner * W
#             T  += chi_inner * W
#         end
#     end

#     norm = X.sumweights^2
#     return real(M1/norm), real(M2/norm), real(Q1/norm), real(Q2/norm), real(Q/norm), real(T/norm), real(norm)
# end

# Rescaled integral from A1,A2 to two normal gaussian with correlation
function inner_trace_legendre(Model::v2_Corrupted_2S_1RSB, X::TI) where TI <: IntegrationMethod
    @extract Model : Oconj params

    J = Oconj[:t] - Oconj[:q] + params[:γ]/2

    # handling zero values of conj params
    hat_q1_clamp = clamp(Oconj[:q1] , 1e-10, Inf) 
    hat_q2_clamp = clamp(Oconj[:q2] , 1e-10, Inf) 
    
    denn = clamp(Oconj[:q1] * Oconj[:q2], 1e-10, Inf)
    ρ = Oconj[:q] / sqrt(denn)

    one_minus_ρ_square = clamp(1 - ρ^2, 1e-10, Inf)


    xi_vals = (-1, 1)

    M1 = M2 = Q1 = Q2 = Q = T = 0.0

    for (i1,z1) in enumerate(X.points)
        w1 = X.weights[i1]
        for (i2,z2) in enumerate(X.points)
            w2 = X.weights[i2]
            g1 = z1
            g2 = ρ*z1 + sqrt(one_minus_ρ_square)*z2

            h1 = Oconj[:m1] + sqrt(hat_q1_clamp) * g1
            h2 = Oconj[:m2] + sqrt(hat_q2_clamp) * g2

            Z = 0.0
            m1i = m2i = chii = 0.0
            for xi1 in xi_vals, xi2 in xi_vals
                w = exp(xi1*h1 + xi2*h2 + J*xi1*xi2)
                Z += w
                m1i += xi1*w
                m2i += xi2*w
                chii += xi1*xi2*w
            end

            m1i /= Z
            m2i /= Z
            chii /= Z

            W = w1 * w2
            M1 += m1i * W
            M2 += m2i * W
            Q1 += m1i^2 * W
            Q2 += m2i^2 * W
            Q  += m1i * m2i * W
            T  += chii * W
        end
    end

    norm = X.sumweights^2

    # global factor 1/2 multiplying each order parameter is missing
    return M1/norm,
           M2/norm,
           Q1/norm,
           Q2/norm,
           Q/norm,
           T/norm,
           norm
end




function lhs!(Model::v2_Corrupted_2S_1RSB, X::TI, Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : aux Oconj
    exprOconj!(Model)
    
    aux["det_cov"] = Oconj[:q1] * Oconj[:q2] - Oconj[:q]^2

    @assert iszero(Onew)

    # Compute order parameters + normalization
    # vals = inner_trace_manual(Model)
    vals = inner_trace_legendre(Model,X)
    
    # Extract order parameters and normalization
    M1, M2, Q1, Q2, Q, T, norm = vals

    aux["norm"] = norm

    # Tolerance threshold for imaginary parts
    imag_tol = 1e-12

    for (name, val) in zip([:m1, :m2, :q1, :q2, :q, :t], (M1, M2, Q1, Q2, Q, T))
        if abs(imag(val)) > imag_tol
            error("Non-negligible imaginary part in $name: $val")
        end
    end

    # Assign only the real parts to Onew
    Onew[:m1] = real(M1)
    Onew[:m2] = real(M2)
    Onew[:q1] = real(Q1)
    Onew[:q2] = real(Q2)
    Onew[:q] = real(Q)
    # Onew[:q]  = max(real(Q), 0.0)
    Onew[:t]  = real(T)

end






########################################################## FREE ENERGY #########################################################

function free_energy_legendre(Model::v2_Corrupted_2S_1RSB, X::TI) where TI <: IntegrationMethod
    @extract Model : params Oconj O

    # Coupling
    J = Oconj[:t] - Oconj[:q] + params[:γ]/2
    ρ = Oconj[:q] / sqrt(Complex(Oconj[:q1] * Oconj[:q2]))

    xi_vals = (-1,1)

    avg_logZ = 0.0

    for (i1,z1) in enumerate(X.points)
        w1 = X.weights[i1]
        for (i2,z2) in enumerate(X.points)
            w2 = X.weights[i2]

            g1 = z1
            g2 = ρ*z1 + sqrt(Complex(1 - ρ^2))*z2

            h1 = Oconj[:m1] + sqrt(Complex(Oconj[:q1])) * g1
            h2 = Oconj[:m2] + sqrt(Complex(Oconj[:q2])) * g2

            Z = 0.0
            for xi1 in xi_vals, xi2 in xi_vals
                Z += exp(xi1*h1 + xi2*h2 + J*xi1*xi2)
            end

            avg_logZ += log(Z) * w1 * w2
        end
    end

    
    avg_logZ = real(avg_logZ)

    # --- Now compute the logdet term using uppercase combinations ---
    β     = params[:β]
    βT    = params[:βT]
    r1    = params[:r1]
    r2    = params[:r2]
    α     = params[:α]
    γ     = params[:γ]

    m1 = O[:m1]
    m2 = O[:m2]
    q1 = O[:q1]
    q2 = O[:q2]
    q  = O[:q]
    t  = O[:t]

    M1 = 1 - β - (βT * β / (1 - βT)) * (r1*m1)^2
    M2 = 1 - β - (βT * β / (1 - βT)) * (r2*m2)^2
    Q1 = -β*q1 - (βT * β / (1 - βT)) * (r1*m1)^2
    Q2 = -β*q2 - (βT * β / (1 - βT)) * (r2*m2)^2
    T  = -r1*r2*β*t - (βT * β / (1 - βT)) * r1*m1*r2*m2
    Q  = -r1*r2*β*q - (βT * β / (1 - βT)) * r1*m1*r2*m2

    denom = (M1-Q1)*(M2-Q2) - (Q-T)^2
    numerator = M2*Q1 + (M1-2*Q1)*Q2 - 2*Q*(T-Q)
    logdet_Xi = real(numerator/denom + log(denom))

    # --- Hatted quadratic term ---
    hat_m1 = Oconj[:m1]
    hat_m2 = Oconj[:m2]
    hat_q1 = Oconj[:q1]
    hat_q2 = Oconj[:q2]
    hat_q  = Oconj[:q]
    hat_t  = Oconj[:t]

    free_term_quadratic = - hat_t * t +
                          0.5*hat_q1*q1 + 0.5*hat_q2*q2 +
                          hat_q*q - hat_m1*m1 - hat_m2*m2 -
                          0.5*hat_q1 - 0.5*hat_q2

    # --- Assemble free energy ---
    F = avg_logZ - (α/2)*logdet_Xi + free_term_quadratic

    return real(-F/params[:β])
end



function free_energy_dummy(Model::v2_Corrupted_2S_1RSB,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj aux 
    exprOconj!(Model)

    # energy_magnetizations = O[:M1] * Oconj[:M1] + O[:M2] * Oconj[:M2] + O[:T] * Oconj[:T] + Oconj[:C]

    # energy_overlaps =  O[:C] * Oconj[:C] + O[:Q] * Oconj[:Q]
    
    # logdetK =  0.5*params[:α] *( log( 1+ params[:β]*( -2 + 2*O[:C] + aux["E2"]*aux["E1"]*params[:β]  ) ) + params[:β]*( -params[:βT]*( O[:M1]^2 * params[:r1]^2 + O[:M2]^2 * params[:r2]^2 + ( (O[:C]-1)*O[:M1]^2 * params[:r1]^2 + (O[:C]-1)*O[:M2]^2 * params[:r2]^2 +2 *(O[:T]-O[:Q])*O[:M1] * O[:M2] * params[:r1]^2 *params[:r2]^2 )*params[:β] ) -2* aux["ΔT"] * ( O[:C] + ( (O[:C]-1)*O[:C] + O[:Q]*( O[:T]-O[:Q] )* params[:r1]^2 *params[:r2]^2 )*params[:β]  ) ) / ( aux["ΔT"]*(1+ params[:β]*( -2 +2*O[:C]+ aux["E1"]*aux["E2"]*params[:β] )   )  ) )
    
    
    # FE = energy_overlaps - energy_magnetizations + aux["Entropic"] - logdetK

    return  - 5
    

end

function free_energy(Model::v2_Corrupted_2S_1RSB, X::TI) where TI <: IntegrationMethod
    
    # --- MANUALLY TOGGLE HERE ---
    
    # return free_energy_legendre(Model, X)
    return free_energy_dummy(Model, X)
    
    # ----------------------------
end

_isSG(Model::v2_Corrupted_2S_1RSB,X::TI,Oold::NamedVec) where TI <: IntegrationMethod = Model.O[:m1] < 1e-5 && Model.O[:m2]<1e-5 && Model.O[:Q] > 1e-4 
_isSR1(Model::v2_Corrupted_2S_1RSB,X::TI,Oold::NamedVec,tol::Float64) where TI <: IntegrationMethod = Model.O[:m1] > 1e-4 

_isSR2(Model::v2_Corrupted_2S_1RSB,X::TI,Oold::NamedVec) where TI <: IntegrationMethod = Model.O[:M2] > 1e-4 && Model.O[:M1] < 1e-5
_isSR(Model::v2_Corrupted_2S_1RSB,X::TI,Oold::NamedVec) where TI <: IntegrationMethod = Model.O[:M2] < 1e-4 && Model.O[:M1] < 1e-4

_isSG1(Model::v2_Corrupted_2S_1RSB,X::TI,Oold::NamedVec,tol::Float64) where TI <: IntegrationMethod = Model.O[:q1] < 1e-4 


_isSignalR(Model::v2_Corrupted_2S_1RSB,X::TI,Oold::NamedVec,tol::Float64) where TI <: IntegrationMethod = Model.O[:m1] > 1e-4 # && Model.O[:m2] > 1e-4



_isSignalNotR(Model::v2_Corrupted_2S_1RSB,X::TI,Oold::NamedVec,tol::Float64) where TI <: IntegrationMethod = Model.O[:m1] < 1e-3



_isParamagnetic(Model::v2_Corrupted_2S_1RSB,X::TI,Oold::NamedVec) where TI <: IntegrationMethod = Model.O[:M1] < 1e-4 && Model.O[:M2] < 1e-4 
