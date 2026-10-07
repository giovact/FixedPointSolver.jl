
struct HomogeneousSpin1Half
    N::Int
    Magnetization::Vector{Int} # eventually simmetrize configurations here
    dos::Vector{BigInt}
    logdos::Vector{FT}
    # occhio che qui sta cosa da errori grandi per N grande (tra exp(logdos) e dos). Quindi boh
end

function HomogeneousSpin1Half(N::Int)
    @assert N>1
    Mvec = collect(range(-N, stop = N, step = 2))
    DOS = collect([dos(N,M) for M in Mvec] )
    #@assert sum(DOS)==2^N # fino a N = 60 va bene, lo hai già controllato. Poi fa casini per i fattoriali suppongo. 
    HomogeneousSpin1Half(N,Mvec, DOS, collect([logdos(N,M)[1] for M in Mvec] ))
end

HomogeneousSpin1Half(Model::TM) where {TM<:FPModel} = HomogeneousSpin1Half(Model.params[:y])

dos(N,M) = BigCombinatorics.Binomial(N,div(N+M,2))
logdos(N,M) = logabsbinomial(N,div(N+M,2))

function recursive_binomial(N,k)
    N==k && return -logfactorial(k)
    return log(N) + recursive_binomial(N-1,k)
end 

CW_energy(N,M,h,J; c = 0.0) = -0.5*J*M^2 - M * h + 0.5*J*N + c
GSenergy(N,h,J) = - 0.5*J*N^2 - N*abs(h) + 0.5*J*N

conf_energylist(C::HomogeneousSpin1Half, h,J) = [CW_energy(C.N, C.Magnetization[i],h,J) for i=1:length(C.Magnetization)]


function trace_CurieWeiss(C::HomogeneousSpin1Half,h,J;return_otherobs = false)
    if C.N == 2
        return trace_2spin_samefield(h,J)
    end
    # # remember, this compute the trace of exp( J \sum_{i<j} σ_i σ_j + h \sum_i σ_i)
    # the coupling J is NOT rescaled by N automatically, so you have to put manually the division factor. 
    Z = 0.0
    m = 0.0
    χ = 0.0

    m2 = 0.0
    mabs = 0.0
    
    mE_gs = - GSenergy(C.N,h,J)
    additive_constant = mE_gs
    for (iM,M) in enumerate(C.Magnetization)
        b = exp(- CW_energy(C.N,M,h,J; c = additive_constant ) + C.logdos[iM] ) 
        Z += b 
        m += b* M
        χ += b* M ^2  # check
        m2 +=b * M^2
        mabs +=b * abs(M)
    end

    return_otherobs && return exp(log(Z) + additive_constant) , m/(C.N*Z) , χ / (Z*C.N*(C.N-1)) - 1/(C.N-1), m2/(C.N^2*Z), mabs/(C.N * Z)  #check actually it is correct
    return exp(log(Z) + additive_constant) , m/(C.N*Z) , χ / (Z*C.N*(C.N-1)) - 1/(C.N-1)
end


function trace_2spin_samefield(h,J)  
    Z = 4*cosh(h)*cosh(h)*cosh(J) + 4*sinh(h)*sinh(h)*sinh(J)
    m = tanh(h + atanh( tanh(J)*tanh(h) ))
    χ = tanh(J + atanh( tanh(h)^2 ))

    return Z, m, χ
end


function trace_2spin(h1,h2,J)  
    Z = 4*cosh(h1)*cosh(h2)*cosh(J) + 4*sinh(h1)*sinh(h2)*sinh(J)
    m1 = tanh(h1 + atanh( tanh(J)*tanh(h2) ))
    m2 = tanh(h2 + atanh( tanh(J)*tanh(h1) ))
    χ = tanh(J + atanh( tanh(h1)*tanh(h2) ))

    return Z, m1, m2, χ
end

###################################################################################################################################

struct TracedSpin1HalfIntegral{TI <: IntegrationMethod} <: IntegrationMethod
    X::TI
    InnerCW::HomogeneousSpin1Half
    Z::Vector{FT}
    m::Vector{FT}
    chi::Vector{FT}
end

TracedSpin1HalfIntegral(Model::TM,X::TI) where {TM <: FPModel,TI <: IntegrationMethod} = TracedSpin1HalfIntegral(X,HomogeneousSpin1Half(Model.params[:y]), zeros(length(X.points)), zeros(length(X.points)), zeros(length(X.points)))

function clean!(X::TracedSpin1HalfIntegral{TI}) where TI <: IntegrationMethod
    fill!(X.Z, 0.0)    
    fill!(X.m, 0.0)
    fill!(X.chi, 0.0)
end

function inner_trace(Model::TM, X) where {TM<:FPModel}
    @assert hasfield(typeof(Model),:InnerCW)
    clean!(X)
    for (iz, z) in enumerate(X.X.points)
        Z, m, χ = inner_trace(Model,z)
        X.Z[iz] = Z
        X.m[iz] = m
        X.chi[iz] = χ
    end
end


### DEPRECATED
### trace

function trace_CurieWeiss(N,h,J)
    ntot = 0    
    Z = 0.0
    m = 0.0
    χ = 0.0

    mE_gs = - GSenergy(N,h,J)
    for M in range(-N, stop = N, step = 2)
        n = Int(dos(N,M))
        b = exp(- CW_energy(N,M,h,J; c = mE_gs))
        ntot += n
        Z += n* b 
        m += n* b* M
        χ += n* b* M ^2  # check
    end
    #@assert ntot == 2^N #debug
    
    N==1 && return Z, m/(N*Z), 1.0

    return Z, m/(N*Z) , χ / (Z*N*(N-1)) - 1/(N-1)
end


tovec(x, n) = [2*((x>>i)&1)-1 for i=0:n-1]

function trace(N,β,h,J)
    #trace

    m=zeros(N)
    C=zeros(N,N)
	Z=0;
    v = zeros(Int,N)
    for x=0:2^N-1
        v .= tovec(x,N)
        p = exp(β* dot(v, (h .+ 0.5 * J*v) ) )
        m .+= v*p
        C .+= (v*v') *p
        Z += p
    end

    m ./= Z
    C ./= Z
    return m,C
end
