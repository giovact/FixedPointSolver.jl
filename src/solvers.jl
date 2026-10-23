abstract type Solver
end

struct FixedPoint <: Solver
    ρ::FT
end

# each function update should return the norm between Onew and Old after the update
function update!(xold,xnew,K::FixedPoint)
    r=norm(xnew - xold,Inf) ;
        xold .*= K.ρ ;
        xold .+= (1.0 - K.ρ) .* xnew;
    return r
end

############################################## Anderson acceleration ########################################
"""
    Anderson(m::Int = 5; ρ = 0.5, λ = 1e-10, restart = 1e3, project = identity)

Anderson acceleration (type II, Walker & Ni 2011) of the fixed-point iteration x -> g(x) = lhs!(Model, ...), to be used
in place of `FixedPoint`, e.g. `solve_SPequations!(Model, X; K = Anderson(5; ρ = 0.5))`. At iteration k, with residuals
f_i = g(x_i) - x_i of the last m+1 iterates,

    γ = argmin |f_k - ΔF γ|₂ (Tikhonov-regularised),   x_{k+1} = x_k + β f_k - (ΔX + β ΔF) γ,   β = 1 - ρ,

where ΔF, ΔX are the differences of consecutive residuals and iterates. With no history (first step, or after a reset)
this is the damped step of `FixedPoint(ρ)`.
- `m`: history length (m = 0 is FixedPoint(ρ)); `ρ`: damping, same convention as FixedPoint; `λ`: relative regularisation.
- `restart`: the history is dropped when |f_k| exceeds `restart` times the smallest residual in the current history.
- `project`: applied to each new iterate (e.g. `x -> clamp.(x, lo, hi)` to keep order parameters in their domain).
- If g(x_k) is not finite (an extrapolated iterate left the domain of the model), the solver goes back to the last iterate
  with a finite image, takes a damped step from it and drops the history; that iteration returns ε = Inf.
- The returned ε is |g(x_k) - x_k|_∞, as for FixedPoint, so `tol` has the same meaning.

The history lives in the solver object and is reset automatically when a new solve starts (detected because the incoming
iterate differs from the one written at the previous call). Do not share one Anderson object between threads. In
`findcriticalpoint` and `scanning`, `K` and `Kfirstrun` must have the same type: pass `Kfirstrun = Anderson(...)` too.
"""
mutable struct Anderson{TP} <: Solver
    m::Int
    ρ::FT
    λ::FT
    restart::FT
    project::TP
    X::Vector{Vector{FT}}     # last iterates x_i
    G::Vector{Vector{FT}}     # their images g(x_i)
    xlast::Vector{FT}         # iterate written at the previous call
    fbest::FT                 # smallest |f|₂ in the current history
    nresets::Int              # number of restarts / recoveries (diagnostic)
end

Anderson(m::Int = 5; ρ::Real = 0.5, λ::Real = 1e-10, restart::Real = 1e3, project = identity) =
    Anderson(m, FT(ρ), FT(λ), FT(restart), project, Vector{FT}[], Vector{FT}[], FT[], Inf, 0)

Base.show(io::IO, K::Anderson) = print(io, "Anderson(m=", K.m, ", ρ=", K.ρ, ")")

function anderson_reset!(K::Anderson)
    empty!(K.X); empty!(K.G)
    K.fbest = Inf
end

function update!(xold,xnew,K::Anderson)
    x, g = Vector{FT}(xold), Vector{FT}(xnew)
    f = g .- x
    r = norm(f, Inf)

    # a new solve started (or the iterate was changed outside the solver): forget the history
    (length(K.xlast) != length(x) || x != K.xlast) && anderson_reset!(K)

    # g(x) not finite: go back to the last iterate with a finite image and take a damped step from there
    if !all(isfinite, g)
        if isempty(K.X)
            xold .*= K.ρ ; xold .+= (1.0 - K.ρ) .* xnew       # nothing to go back to: behave as FixedPoint
            return r
        end
        xn = K.project(K.ρ .* K.X[end] .+ (1 - K.ρ) .* K.G[end])
        anderson_reset!(K); K.nresets += 1
        xold .= xn; K.xlast = Vector{FT}(xold)
        return Inf
    end

    push!(K.X, x); push!(K.G, g)
    length(K.X) > K.m + 1 && (popfirst!(K.X); popfirst!(K.G))
    fn = norm(f)
    if fn > K.restart * K.fbest                                # residual exploded: restart from the current pair
        K.X, K.G = [x], [g]; K.fbest = Inf; K.nresets += 1
    end
    K.fbest = min(K.fbest, fn)

    mk = length(K.X) - 1
    β = 1 - K.ρ
    if mk == 0
        xn = x .+ β .* f
    else
        F  = [K.G[i] .- K.X[i] for i in 1:mk+1]
        ΔF = reduce(hcat, [F[i+1] .- F[i] for i in 1:mk])
        ΔX = reduce(hcat, [K.X[i+1] .- K.X[i] for i in 1:mk])
        A = ΔF' * ΔF
        γ = (A + (K.λ * tr(A) / mk + floatmin(FT)) * I) \ (ΔF' * f)
        xn = x .+ β .* f .- (ΔX .+ β .* ΔF) * γ
        all(isfinite, xn) || (xn = x .+ β .* f; anderson_reset!(K); K.nresets += 1)
    end
    xn = K.project(xn)

    xold .= xn
    K.xlast = Vector{FT}(xold)
    return r
end
