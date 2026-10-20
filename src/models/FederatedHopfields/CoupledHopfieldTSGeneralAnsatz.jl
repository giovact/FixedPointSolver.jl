
# SDN model with homogeneous corruption r, GENERAL self/cross RS ansatz
# (dAT_SDN_v2 notes, Sec. "General self/cross ansatz: equations to implement"):
#
#   Q_aa = I + t (1 - I),    Q_ab = q 1 + (qt - q) I   (a != b),    m^a_u = m
#
# qt = self-overlap (same student, two replicas), q = cross-overlap (two students, two replicas).
# For qt = q this is the ansatz (B11) of Corrupted.jl; at r = 1 the two coincide.
#
# The local block has a common field and one private field per student,
#   h_u = mhat + sqrt(qhat) z + sqrt(qthat - qhat) z_u,     J = β^ν γ/y + that - qhat,
# so every average is over z (the outer quadrature X.X of the library) and over the y private
# Gaussians z_u (inner quadrature of PrivateNoiseIntegral below).

struct CoupledHopfieldTSGeneralAnsatz <: FPModel
    params::NamedTuple{(:β,:α,:r,:βT,:γ,:y,:ν),Tuple{FT,FT,FT,FT,FT,Int, Bool}}
    O::NamedVec      # order parameters  -> M, P, Q (cross q), Qt (self qt), T
    Oconj::NamedVec  # conjugates        -> M, Q (cross qhat), Qt (self qthat), T
    InnerCW::HomogeneousSpin1Half
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors
CoupledHopfieldTSGeneralAnsatz() = CoupledHopfieldTSGeneralAnsatz(.1, .1, 1.0, .1, 1.0, 2, false, ones(5)) # dummy constructor for testing
CoupledHopfieldTSGeneralAnsatz(β,α,r,βT,γ,y,ν) = CoupledHopfieldTSGeneralAnsatz(β,α,r,βT,γ,y,ν,zeros(5))
CoupledHopfieldTSGeneralAnsatz(β,α,r,βT,γ,y,ν,O) = CoupledHopfieldTSGeneralAnsatz((β=β,α=α,r=r,βT=βT,γ=γ,y=y,ν=ν),NamedArray(O,[:M;:P;:Q;:Qt;:T]),NamedArray(zeros(4),[:M;:Q;:Qt;:T]), HomogeneousSpin1Half(y),Dict{String,FT}())

create_model(Model::CoupledHopfieldTSGeneralAnsatz, params) = CoupledHopfieldTSGeneralAnsatz(params,Model.O, Model.Oconj, HomogeneousSpin1Half(Int(params[:y])),Dict{String,Float64}())


function sanity_checks!(Model::CoupledHopfieldTSGeneralAnsatz)
    Model.params[:βT] >= 1 && throw("You should use a βT < 1 !!!!")
    Model.params[:y] < 2 && throw("Use y >= 2 (y = 1 is HopfieldTS.jl)")
end

########################################################## PRIVATE-NOISE INTEGRATION #########################################################
# normalised weights of a 1-d quadrature of the library (HermiteQuadrature, LegendreQuadrature, UniformInterp):
# w_i = integrate(Xp, e_i), so that sum_i w_i f(x_i) = integrate(Xp, f)
quadrature_weights(Xp::TI) where TI <: IntegrationMethod = [integrate(Xp, FT.(1:length(Xp.points) .== i)) for i in 1:length(Xp.points)]

# All the statistics used below are symmetric under permutations of the students, so the tensor
# grid of the y private Gaussians can be reduced to its multisets (non-decreasing node tuples),
# each weighted by multinomial(y; counts) * prod(w): C(n+y-1, y) points instead of n^y.
function private_nodes(y::Int, Xp::TI) where TI <: IntegrationMethod
    x, w = Xp.points, quadrature_weights(Xp)
    n = length(x)
    cols = Vector{Vector{FT}}(); wts = FT[]
    idx = ones(Int, y)
    while true
        counts = zeros(Int, n); for i in idx; counts[i] += 1; end
        push!(cols, x[idx])
        push!(wts, exp(logfactorial(y) - sum(logfactorial.(counts))) * prod(w[idx]))
        k = y
        while k >= 1 && idx[k] == n; k -= 1; end
        k == 0 && break
        idx[k] += 1
        for j in k+1:y; idx[j] = idx[k]; end
    end
    return reduce(hcat, cols), wts ./ sum(wts)
end

# default number of private nodes per dimension: largest n <= 40 with C(n+y-1, y) <= 4000
function default_npriv(y::Int)
    n = 3
    while n < 40 && binomial(n + y, y) <= 4000; n += 1; end
    return n
end

struct PrivateNoiseIntegral{TI <: IntegrationMethod, TP <: IntegrationMethod} <: IntegrationMethod
    X::TI                  # quadrature over the common Gaussian z (X.points)
    Xpriv::TP              # 1-d quadrature used for each private Gaussian z_u
    configs::Matrix{FT}    # 2^y x y spin configurations
    S2::Vector{FT}         # (sum_u ξ_u)^2 of each configuration
    zp::Matrix{FT}         # y x K private nodes
    wp::Vector{FT}         # K weights, sum to 1
    # private-noise averages at each outer point z (length(X.points))
    logZ::Vector{FT}       # < ln Z >
    m::Vector{FT}          # < M_u >
    qt::Vector{FT}         # < M_u^2 >
    q::Vector{FT}          # < M_u M_v >,  u != v
    t::Vector{FT}          # < C_uv >,     u != v  (full, not connected)
    # dAT statistics (filled only by inner_trace(...; dat = true))
    Gpp::Vector{FT}; Gpq::Vector{FT}; Gqq::Vector{FT}
    A11::Vector{FT}; A12::Vector{FT}; A22::Vector{FT}; Ysum::Vector{FT}
end

# Xpriv: any 1-d Gaussian quadrature of the library with `points` (default: HermiteQuadrature)
function PrivateNoiseIntegral(Model::CoupledHopfieldTSGeneralAnsatz, X::TI;
                              Xpriv::TP = HermiteQuadrature(default_npriv(Model.params[:y]))) where {TI <: IntegrationMethod, TP <: IntegrationMethod}
    y = Model.params[:y]
    configs = FT[2 * ((c >> (u-1)) & 1) - 1 for c in 0:2^y-1, u in 1:y]
    zp, wp = private_nodes(y, Xpriv)
    np = length(X.points)
    PrivateNoiseIntegral(X, Xpriv, configs, vec(sum(configs, dims = 2)) .^ 2, zp, wp, [zeros(np) for _ in 1:12]...)
end

print_integration_params(X::PrivateNoiseIntegral) = string(print_integration_params(X.X), "_priv", print_integration_params(X.Xpriv))

# wraps the library quadrature X (over the common z) with the private-noise quadrature; pass an already built
# PrivateNoiseIntegral(Model, X; Xpriv = HermiteQuadrature(n)) to choose the private nodes yourself
# (if y has changed in the meantime, e.g. after create_model, the private grid is rebuilt with the default size)
set_integrationmethod(Model::CoupledHopfieldTSGeneralAnsatz,X::TI) where TI <: IntegrationMethod = PrivateNoiseIntegral(Model,X)
set_integrationmethod(Model::CoupledHopfieldTSGeneralAnsatz,X::PrivateNoiseIntegral) = size(X.configs, 2) == Model.params[:y] ? X : PrivateNoiseIntegral(Model, X.X)

########################################################## SELF-CONSISTENT EQUATIONS #########################################################
DeltaT(Model::CoupledHopfieldTSGeneralAnsatz) = 1 - Model.params[:βT]
# eigenvalues of D = I - β R∘(Q_aa - Q_ab):  Δ (on Π_⊥) and Δy (on Π_∥)
Delta(Model::CoupledHopfieldTSGeneralAnsatz)  = 1 - Model.params[:β] + Model.params[:β] * Model.O[:Qt] + Model.params[:β] * Model.params[:r]^2 * (Model.O[:T] - Model.O[:Q])
Deltay(Model::CoupledHopfieldTSGeneralAnsatz) = 1 - Model.params[:β] + Model.params[:β] * Model.O[:Qt] - Model.params[:β] * Model.params[:r]^2 * (Model.O[:T] - Model.O[:Q]) * (Model.params[:y] - 1)

is_free_energy_definite(Model::CoupledHopfieldTSGeneralAnsatz) = (Delta(Model)>0) && (Deltay(Model)>0)

function Deltas!(Model::CoupledHopfieldTSGeneralAnsatz)
    @extract Model : params O aux
    β, r, y = params[:β], params[:r], params[:y]
    aux["ΔT"] = DeltaT(Model)
    aux["Δ"]  = clamp(Delta(Model), 1e-30, Inf)
    aux["Δy"] = clamp(Deltay(Model), 1e-30, Inf)
    # eigenvalues of F = β R∘Q_ab + (βT β r^2 m^2/ΔT) 1
    aux["f_par"]  = β * (O[:Qt] + (y-1) * r^2 * O[:Q]) + y * params[:βT] * β * r^2 * O[:M]^2 / aux["ΔT"]
    aux["f_perp"] = β * (O[:Qt] - r^2 * O[:Q])
end

#### Conjugate order params SP (stationarity of the energetic term)
function exprOconj!(Model::CoupledHopfieldTSGeneralAnsatz)
    @extract Model : O Oconj params aux
    Deltas!(Model)
    β, βT, r, α, y = params[:β], params[:βT], params[:r], params[:α], params[:y]
    Δ, Δy, ΔT = aux["Δ"], aux["Δy"], aux["ΔT"]
    fp, fo = aux["f_par"], aux["f_perp"]

    Oconj[:M]  = α * βT * β * r^2 * O[:M] / (ΔT * Δy)
    Oconj[:Qt] = (α * β / y) * ( fp / Δy^2 + (y-1) * fo / Δ^2 )          # self
    Oconj[:Q]  = (α * β * r^2 / y) * ( fp / Δy^2 - fo / Δ^2 )             # cross
    Oconj[:T]  = (α * β * r^2 / y) * ( 1 / Δy - 1 / Δ ) + Oconj[:Q]       # that - qhat = α β^2 r^4 (t-q)/(Δ Δy)
end

# fields of the local block
Jeff(Model::CoupledHopfieldTSGeneralAnsatz) = (Model.params[:β]^Model.params[:ν]) * Model.params[:γ] / Model.params[:y] + (Model.Oconj[:T] - Model.Oconj[:Q])
sqrt_qhat(Model::CoupledHopfieldTSGeneralAnsatz) = sqrt(max(Model.Oconj[:Q], 0.0))                    # common field
sqrt_phat(Model::CoupledHopfieldTSGeneralAnsatz) = sqrt(max(Model.Oconj[:Qt] - Model.Oconj[:Q], 0.0)) # private fields

# local y-spin block with heterogeneous fields h, coupling J: fills p (Boltzmann weights) and M, returns ln Z
function local_block!(p, M, configs, S2, h, J, y)
    mul!(p, configs, h)
    p .+= 0.5 * J .* (S2 .- y)
    c = maximum(p)
    p .= exp.(p .- c)
    Z = sum(p)
    p ./= Z
    mul!(M, transpose(configs), p)
    return log(Z) + c
end

function inner_trace(Model::CoupledHopfieldTSGeneralAnsatz, X::PrivateNoiseIntegral; dat::Bool = false)
    @extract Model : params Oconj
    y = params[:y]
    J = Jeff(Model)
    a, b = sqrt_qhat(Model), sqrt_phat(Model)
    Model.aux["qhat_negative"] = Oconj[:Q] < 0     # should never be true; the common field needs qhat >= 0
    K = b > 0 ? length(X.wp) : 1                    # no private noise -> a single private node
    nc = size(X.configs, 1)
    p, M, h = zeros(nc), zeros(y), zeros(y)
    C, χ, Pq, χb, B = zeros(y, y), zeros(y, y), Matrix{FT}(I, y, y) .- 1 / y, zeros(y, y), zeros(y, y)
    for (iz, z) in enumerate(X.X.points)
        lz = mm = m2 = mc = tc = 0.0
        gpp = gpq = gqq = a11 = a12 = a22 = ys = 0.0
        for k in 1:K
            wk = K == 1 ? 1.0 : X.wp[k]
            for u in 1:y; h[u] = Oconj[:M] + a * z + (K == 1 ? 0.0 : b * X.zp[u, k]); end
            lz += wk * local_block!(p, M, X.configs, X.S2, h, J, y)
            sM, sM2 = sum(M), sum(abs2, M)
            mm += wk * sM / y
            m2 += wk * sM2 / y
            mc += wk * (sM^2 - sM2) / (y * (y - 1))
            tc += wk * (dot(p, X.S2) - y) / (y * (y - 1))
            if dat
                # connected correlation χ = C - M M^T of this sample, and the sector invariants of the notes
                mul!(C, transpose(X.configs), X.configs .* p)
                χ .= C .- M * M'
                v = vec(sum(χ, dims = 2)); s = sum(v); v2 = sum(abs2, v)
                τ, τ2 = tr(χ), sum(abs2, χ)
                gpp += wk * s^2 / y^2
                gpq += wk * (v2 / y - s^2 / y^2)
                gqq += wk * (τ2 - 2 * v2 / y + s^2 / y^2)
                av = v .- s / y                                   # Π_⊥ χ 1
                a11 += wk * (2 * sum(abs2, av) + 2 * s * (τ - s / y))
                mul!(B, Pq, χ); mul!(χb, B, Pq)                   # B = Π_⊥ χ, χb = Π_⊥ χ Π_⊥
                a22 += wk * (sum(abs2, diag(χb)) - sum(abs2, χb) / y)
                a12 += wk * 2 * sum(av[i] * B[i, j] * Pq[i, j] for i in 1:y, j in 1:y)
                ys  += wk * 0.5 * (τ^2 + τ2)                      # trace of V -> χVχ on Sym(y)
            end
        end
        X.logZ[iz], X.m[iz], X.qt[iz], X.q[iz], X.t[iz] = lz, mm, m2, mc, tc
        if dat
            X.Gpp[iz], X.Gpq[iz], X.Gqq[iz] = gpp, gpq, gqq
            X.A11[iz], X.A12[iz], X.A22[iz], X.Ysum[iz] = a11, a12, a22, ys
        end
    end
end

function lhs!(Model::CoupledHopfieldTSGeneralAnsatz,X::PrivateNoiseIntegral,Onew::NamedVec)
    exprOconj!(Model)
    @assert iszero(Onew)
    inner_trace(Model, X)
    Onew[:M]  = integrate(X.X, X.m)
    Onew[:P]  = 0.0
    Onew[:Q]  = integrate(X.X, X.q)
    Onew[:Qt] = integrate(X.X, X.qt)
    Onew[:T]  = integrate(X.X, X.t)
end

########################################################## FREE ENERGY #########################################################
# -β f_RS = Φ / y, with
# Φ = - y(y-1)/2 that t + y/2 qthat (qt - 1) + y(y-1)/2 qhat q - y mhat m
#     - α/2 [ ln Δy + (y-1) ln Δ - f_par/Δy - (y-1) f_perp/Δ ] + < ln Z >_{z, z_u}
# at the current O and Oconj (no recomputation of the conjugates; used for the stationarity checks)
function Phi_RS(Model::CoupledHopfieldTSGeneralAnsatz, X::PrivateNoiseIntegral)
    @extract Model : params O Oconj aux
    Deltas!(Model)
    α, y = params[:α], params[:y]
    inner_trace(Model, X)
    overlaps = - 0.5 * y * (y-1) * Oconj[:T] * O[:T] + 0.5 * y * Oconj[:Qt] * (O[:Qt] - 1) +
               0.5 * y * (y-1) * Oconj[:Q] * O[:Q] - y * Oconj[:M] * O[:M]
    logdet = log(aux["Δy"]) + (y-1) * log(aux["Δ"]) - aux["f_par"] / aux["Δy"] - (y-1) * aux["f_perp"] / aux["Δ"]
    entropy = integrate(X.X, X.logZ)
    aux["Phi_overlaps"], aux["Phi_logdet"], aux["Phi_entropy"] = overlaps, -0.5 * α * logdet, entropy
    return overlaps - 0.5 * α * logdet + entropy
end

# free energy per student, f = -Φ/(β y)
function free_energy(Model::CoupledHopfieldTSGeneralAnsatz,X::TI) where TI <: IntegrationMethod
    exprOconj!(Model)
    return - Phi_RS(Model, X) / (Model.params[:β] * Model.params[:y])
end

# same conventions as Corrupted.jl; the spin glass is detected on the SELF-overlap Qt: in the general ansatz a
# spin glass can have Q = 0 < Qt (e.g. at γ = 0, where the students freeze independently)
_isSG(Model::CoupledHopfieldTSGeneralAnsatz,X::TI,Oold::NamedVec,tol::Float64) where TI <: IntegrationMethod = Model.O[:Qt] < 1e-5 && Model.O[:P] < 1e-5
_isSignalR(Model::CoupledHopfieldTSGeneralAnsatz,X::TI,Oold::NamedVec,tol::Float64) where TI <: IntegrationMethod = Model.O[:M] > 1e-2
_isExampleR(Model::CoupledHopfieldTSGeneralAnsatz,X::TI,Oold::NamedVec,tol::Float64) where TI <: IntegrationMethod = Model.O[:P] > 1e-4 && Model.O[:M] < 1e-5
_isParamagnetic(Model::CoupledHopfieldTSGeneralAnsatz,X::TI,Oold::NamedVec,tol::Float64) where TI <: IntegrationMethod = Model.O[:M] < 1e-3

########################################################## RS (dAT) STABILITY #########################################################
"""
    eigs_dat_general!(Model::CoupledHopfieldTSGeneralAnsatz, X::PrivateNoiseIntegral)

Replicon stability of the general (self/cross) RS solution (dAT_SDN_v2 notes, section
"General self/cross ansatz: equations to implement", subsection "Replicon (dAT) instability"). Stored in Model.aux:
- "Gs_11", "Gs_12", "Gs_21", "Gs_22": entropic matrix Γ^s (diagonal for the ansatz (B11))
- "Lambda_S", "Lambda_A", "Lambda_Y": largest eigenvalue of each sector (Lambda_Y = -Inf if y < 4)
- "dAT_S_det" = (Δy Δ)^2 (1 - Tr Ts + det Ts), "dAT_S_tr" = (Δy Δ)^2 (2 - Tr Ts): division-free margins
- "max_Lambda" = max of the three, "RS_valid", "is_RS_stable": stable iff all Λ < 1 (equivalently all margins > 0)
Called by compute_aux! at convergence of solve_SPequations!; NaN-safe (Λ = NaN if the solution is not finite).
"""
# largest real eigenvalue, NaN if the matrix is not finite (eigvals would throw)
maxeig(A) = all(isfinite, A) ? maximum(real.(eigvals(A))) : NaN

function eigs_dat_general!(Model::CoupledHopfieldTSGeneralAnsatz, X::PrivateNoiseIntegral)
    @extract Model : aux params
    exprOconj!(Model)
    inner_trace(Model, X; dat = true)
    α, β, r, y = params[:α], params[:β], params[:r], params[:y]
    ab2 = α * β^2
    Δ, Δy = Delta(Model), Deltay(Model)          # unclamped
    aux["RS_valid"] = (Δ > 0) && (Δy > 0) && all(isfinite, X.m)

    cp, co, g = 1 + (y-1) * r^2, 1 - r^2, y - 1 + r^2

    # (S) symmetric sector: Ts = αβ^2 Γs Hs diag(Δy^-2, Δ^-2) Hs
    Gpp, Gpq, Gqq = integrate(X.X, X.Gpp), integrate(X.X, X.Gpq), integrate(X.X, X.Gqq)
    Γs = [Gpp Gpq; Gpq/(y-1) Gqq/(y-1)]
    Hs = [cp (y-1)*co; co g] ./ y
    Ts = ab2 .* (Γs * Hs * Diagonal([1/Δy^2, 1/Δ^2]) * Hs)
    aux["Gs_11"], aux["Gs_12"], aux["Gs_21"], aux["Gs_22"] = Γs[1,1], Γs[1,2], Γs[2,1], Γs[2,2]
    aux["Lambda_S"] = maxeig(Ts)
    # division-free margins: multiply by (Δy Δ)^2 (det Hs = r^2)
    Ms = Hs * Diagonal([Δ^2, Δy^2]) * Hs                  # = (Δy Δ)^2 Hs diag(Δy^-2, Δ^-2) Hs
    trS = ab2 * tr(Γs * Ms)
    aux["dAT_S_det"] = (Δy * Δ)^2 - trS + ab2^2 * (Γs[1,1] * Γs[2,2] - Γs[1,2] * Γs[2,1]) * r^4
    aux["dAT_S_tr"]  = 2 * (Δy * Δ)^2 - trS

    # (A) antisymmetric sector
    A11, A12, A22 = integrate(X.X, X.A11), integrate(X.X, X.A12), integrate(X.X, X.A22)
    if y == 2
        aux["Lambda_A"] = ab2 * (A11 / (2 * y * (y-1))) / (Δy * Δ)
    else
        c1, c2 = 2 * y, (y - 2) / y
        Γa = [A11/((y-1)*c1) A12/((y-1)*c1); A12/((y-1)*c2) A22/((y-1)*c2)]
        Ha = [(2 + (y-2)*r^2)/y (y-2)*co/y^2; 2*co (y-2+2*r^2)/y]
        Ta = ab2 .* (Γa * Ha * Diagonal([1/(Δy*Δ), 1/Δ^2]) * Ha)
        aux["Lambda_A"] = maxeig(Ta)
    end

    # (Y) (y-2,2) sector, y >= 4: Γ_Y = [Tr_Sym - Tr_S - Tr_A] / dim
    if y >= 4
        ΓY = (integrate(X.X, X.Ysum) - Gpp - Gqq/(y-1) - A11/(2*y) - A22*y/(y-2)) / (y*(y-3)/2)
        aux["Lambda_Y"] = ab2 * r^4 * ΓY / Δ^2
    else
        aux["Lambda_Y"] = -Inf
    end

    aux["max_Lambda"] = max(aux["Lambda_S"], aux["Lambda_A"], aux["Lambda_Y"])    # NaN if any is NaN
    aux["is_RS_stable"] = (aux["RS_valid"] == 1.0) && aux["max_Lambda"] < 1.0
    return nothing
end

_is_RS_stable(Model::CoupledHopfieldTSGeneralAnsatz,X::TI,Oold::NamedVec,tol::Float64) where TI <: IntegrationMethod = Model.aux["max_Lambda"] < 1.0

# called by solve_SPequations! at convergence (after free_energy!), as for Corrupted.jl
function compute_aux!(Model::CoupledHopfieldTSGeneralAnsatz,X::PrivateNoiseIntegral)
    eigs_dat_general!(Model, X)
    SG_signal_instability!(Model)
end

# Instability of an SG fixed point (M = 0) towards the signal (single student: HopfieldTS.jl; notes dAT_SDN/HopfieldTS_1RSB, Sec. 7).
# A perturbation δm_u creates the fields δm̂_u = (α βT β r²/ΔT) Σ_v (D⁻¹)_uv δm_v, D = Δy Π_∥ + Δ Π_⊥, and responds through the
# averaged connected correlations of the local block, <χ> = χ_∥ Π_∥ + χ_⊥ Π_⊥ (average over z and the private fields z_u):
#   λm_sym = (α βT β r²/ΔT) χ_∥/Δy,   λm_anti = (α βT β r²/ΔT) χ_⊥/Δ,
# with, by the sum rules of the block, χ_∥ = 1 + (y-1)T - Qt - (y-1)Q and χ_⊥ = 1 - Qt - T + Q. λm > 1: M = 0 unstable
# (only meaningful at M = 0). At γ = 0 the students are independent and λm_sym is the single-student loop gain.
function SG_signal_instability!(Model::CoupledHopfieldTSGeneralAnsatz)
    @extract Model : params O aux
    y = params[:y]
    pre = params[:α] * params[:βT] * params[:β] * params[:r]^2 / DeltaT(Model)
    aux["χ_∥"] = 1 + (y - 1) * O[:T] - O[:Qt] - (y - 1) * O[:Q]
    aux["χ_⊥"] = 1 - O[:Qt] - O[:T] + O[:Q]
    aux["λm_sym"] = pre * aux["χ_∥"] / Deltay(Model)
    aux["λm_anti"] = pre * aux["χ_⊥"] / Delta(Model)
    aux["λm"] = max(aux["λm_sym"], aux["λm_anti"])
    aux["is_SG_unstable"] = aux["λm"] > 1
end

# true while the fixed point is a stable SG (keepgoing function for critical_line.jl: scanning α at fixed β stops at the SG -> signal line)
_isSGstable(Model::CoupledHopfieldTSGeneralAnsatz,X::TI,Oold::NamedVec,tol::Float64) where TI <: IntegrationMethod = Model.O[:M] < 1e-3 && Model.aux["λm"] < 1
