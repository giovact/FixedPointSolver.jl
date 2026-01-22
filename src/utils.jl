gofy(y, resc) = resc ? y : 1
tovec(x, n) = [2*((x>>i)&1)-1 for i=0:n-1]

Herf(x) = 0.5 * erfc(sqrt(0.5) * x)
logHerf(x) = logerfc(sqrt(0.5)*x) - log(2.0)
log2cosh(x) = abs(x) + log1p(exp(-2abs(x)))
logcosh(x) = log2cosh(x) - log(2.0)

threshold_to_0(x,ε) = ifelse(abs(x)<ε,0,x)

"""
G(x) = exp(-x^2/2) / F(√(2π))
H(x) = erfc(x /F(√2)) / 2
GH(x) = 2 / erfcx(x/F(√2)) / F(√(2π))
HG(x) =  F(√(2π))*erfcx(x/F(√2)) / 2
G(x, Δ) = G(x, 0, Δ)
G(x, μ, Δ) = exp(-(x-μ)^2/(2Δ)) / √(2π*Δ)
logH(x) = sf_log_erfc(x/√2) - log(2)
logG(x, μ, Δ) = -(x-μ)^2/(2Δ) - log(2π*Δ)/2
logG(x) = -x^2/2 - log(2π)/2

lrelu(x, γ=0.1f0) = max(x, γ*x)
log2cosh(x) = abs(x) + log1p(exp(-2abs(x)))
logcosh(x) = log2cosh(x) - log(2)

θfun(x) = x > 0 ? 1 : 0
"""

compute_aux!(Model::TM,X::TI;setNaNs::Bool = false) where {TM <: FPModel, TI <: IntegrationMethod} = return

function force_op!(O::NamedVec,force_dict::Dict)
    for op in keys(force_dict)
        O[op] = force_dict[op]
    end
end

function force_op!(Model::TM,force_dict::Dict ) where {TM <: FPModel}
    for op in keys(force_dict)
        Model.O[op] = force_dict[op]
    end
end

initialize_history(Model::TM) where TM <: FPModel = MVHistory()

function update_history!(Model::TM, history::TH) where {TM <: FPModel, TH <:MVHistory}
    for n in names(Model.O,1)
        push!(history,n,Model.O[n])
    end
    if hasfield(typeof(Model), :Oconj)
        for n in names(Model.Oconj,1)
            push!(history,Symbol(n,"conj"),Model.Oconj[n])
        end
    end
end

function reset_orderparams!(O,tol)
    for (i,x) in enumerate(O)
        if O[i]!=0.0
            O[i] = max(O[i],tol)
        end
    end
end

# returns a NamedTuple equal to params but with the key k changed to value v
reset_parameter(params,k,v) = @assert k in keys(params) &&  return (; params..., k=>v)

function params_to_textstring(a::Dict)
    out = ""
    for k in keys(a)
        out *= string(k,"_",a[k], "_")
    end
    return out
end

function params_to_textstring(Model::TM;separator::String = "_",exclude::Vector{Symbol}= Symbol[],equalsymbol::String=separator) where TM <: FPModel
    out = ""; 
    for k in keys(Model.params)
        if k ∉ exclude
            out *= string(k,equalsymbol,Model.params[k], separator)
        end
    end
    out
end

function params_to_textstring(params::NamedTuple;separator::String = "_",exclude::Vector{Symbol}= Symbol[],equalsymbol::String=separator)
    out = ""; 
    for k in keys(Model.params)
        if k ∉ exclude
            out *= string(k,equalsymbol,params[k], separator)
        end
    end
    out
end


function print_params(Model::TM;separator::String = "_",exclude::Vector{Symbol}= Symbol[]) where TM <: FPModel 
    string( typeof(Model), " - ", params_to_textstring(Model;separator = separator, exclude = exclude))
end


function print_tuple_params(Model::TM;exclude::Vector{Symbol}= Symbol[]) where TM <: FPModel 
    vec_params = []
    for k in keys(Model.params)
        if k ∉ exclude
            push!(vec_params,(;k=>Model.params[k]))
        end
    end
    return merge(vec_params...)
end


##### greek letter mapping #####
greek_to_latin = Dict(
        "α" => "alpha","β" => "beta", "γ" => "gamma", "δ" => "delta", "ε" => "epsilon", "ζ" => "zeta", "η" => "eta", "θ" => "theta","ι" => "iota",
        "κ" => "kappa", "λ" => "lambda", "μ" => "mu", "ν" => "nu", "ξ" => "xi", "ο" => "omicron", "π" => "pi", "ρ" => "rho", "σ" => "sigma",
        "τ" => "tau", "υ" => "upsilon", "φ" => "phi", "χ" => "chi", "ψ" => "psi","ω" => "omega",
        "Α" => "Alpha", "Β" => "Beta", "Γ" => "Gamma", "Δ" => "Delta", "Ε" => "Epsilon", "Ζ" => "Zeta", "Η" => "Eta", "Θ" => "Theta", "Ι" => "Iota", 
        "Κ" => "Kappa", "Λ" => "Lambda", "Μ" => "Mu", "Ν" => "Nu", "Ξ" => "Xi", "Ο" => "Omicron", "Π" => "Pi", "Ρ" => "Rho", "Σ" => "Sigma", "Τ" => "Tau", 
        "Υ" => "Upsilon","Φ" => "Phi", "Χ" => "Chi", "Ψ" => "Psi", "Ω" => "Omega"
)

latin_to_greek= Dict(value => key for (key, value) in greek_to_latin)

function replace_Greek(in::String; mapping = greek_to_latin)

    for greekletter in keys(mapping)
        if occursin(greekletter,in) 
            in = replace(in, greekletter => mapping[greekletter] )
        end
    end
    return in
end

dumb(T::DataType) = T <: Integer ? 2 : zero(T)

function collect_functions_with_underscore(mod::Module)
    func_dict = Dict()
    for name in names(mod, all=true)
        str_name = string(name)
        if startswith(str_name, "_") && isdefined(mod, name) && isa(getfield(mod, name), Function)
            # Strip the underscore and map the name to the function object
            key = replace(str_name, "_" => "", count=1)
            func_dict[key] = getfield(mod, name)
        end
    end
    return func_dict
end



function create_subdir(dir::String,subdir::String)
    !ispath(dir) && mkpath(dir)
    !ispath(string(dir,"/",subdir)) && (subdir!="") && mkpath(string(dir,"/",subdir))
end

suffix(in::String) = ifelse(isempty(in), "",string(in,"_"))