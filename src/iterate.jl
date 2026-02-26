struct New
    O::NamedVec 
    Oconj::NamedVec
end

function New(Model::TM) where {TM<: FPModel}
    O = similar(Model.O)
    Oconj = n_conj_order_params(Model)==0 ? NamedArray([0.0],[:empty]) : similar(Model.Oconj)
    return New(O,Oconj)
end

function clear!(a::New)
    fill!(a.O, 0.0)
    fill!(a.Oconj, 0.0)
end

function update!(old,new,ρ)
    r=norm(new-old,Inf);
        old .*=ρ ;
        old .+= (1-ρ).*new;
    return r
end

function update_SPequations!(Model::TM,Onew::NamedVec,damp::FT) where TM <:FPModel
    return update!(Model.O,Onew,damp)
end

function wrap_output(flag, vareps, iter, dict_hist, savehist)
    savehist && return flag, iter, vareps, dict_hist
    return flag, iter, vareps
end

enforce_initial_condition!(Model::TM)  where {TM <: FPModel} = nothing
sanity_checks!(Model::TM) where {TM <: FPModel} = nothing

function solve_SPequations!(Model::TM,Xstart::TI;
                            niter::Int = 100000,
                            K::TU = FixedPoint(0.9),
                            tol::FT = 1e-4,
                            savehistory::Bool = false,
                            showinfo::Bool=true,
                            printprogress::Bool=true,
                            printaux::Vector{String} = String[],
                            force_dict = Dict()) where {TM <: FPModel, TI <: IntegrationMethod,  TU<:Solver}
    
    

    sanity_checks!(Model)
    force_op!(Model,force_dict)
    enforce_initial_condition!(Model)
    history = initialize_history(Model)
    X = set_integrationmethod(Model,Xstart)

    new_OP = New(Model)
    ε = 10.0
    prog = printprogress ? ProgressUnknown(string("FP equations on ", print_params(Model; separator= " - "), " - iter = ")) : nothing
    for iter=1:niter
        clear!(new_OP)
        lhs!(Model,X,new_OP.O)
        force_op!(new_OP.O,force_dict) # così se lo calcola in lhs ma poi lo sovrascrive. va bene comunque per ora
        ε = update!(Model.O,new_OP.O,K)
        printprogress && ProgressMeter.next!(prog; showvalues = printvalues(Model,ε, printaux))
        savehistory && update_history!(Model,history)
        if ε < tol
            compute_stuff!(Model,X)
            showinfo && @info string(" Converged in $iter iterations -> error = $ε" )
            printprogress && ProgressMeter.finish!(prog)
            return wrap_output(1,ε, iter, history, savehistory)
        end                                                                                                                                                                                            
        # faster end if NaN is encountered
        if any(isnan.(Model.O))
            showinfo && @info string(" found NaNs at iter $iter")
            compute_stuff!(Model,X;setNaNs = true)
            printprogress && ProgressMeter.finish!(prog)
            return wrap_output(-1,ε, iter, history, savehistory)
        end
    end
    compute_stuff!(Model,X)
    printprogress && ProgressMeter.finish!(prog)
    showinfo && @info string(" NOT Converged in $niter iterations with $K . Final ε= $ε ")

    return wrap_output(0,ε, niter, history, savehistory)
end

function compute_stuff!(Model::TM,X::TI;setNaNs::Bool = false) where {TM<: FPModel, TI <: IntegrationMethod}
    if setNaNs 
        Model.aux["f"] = NaN
        compute_aux!(Model,X)
    else
        free_energy!(Model,X)
        compute_aux!(Model,X)
    end
end

function free_energy!(Model::TM,X::TI) where {TM<: FPModel, TI <: IntegrationMethod}
    Model.aux["f"] = free_energy(Model,X)
end



function solve_Multiple_initialconditions(Model::TM, X::TI,O0vec::Vector{Vector{FT}}; # generalize to multiple Op but ok 
                                                niter::Int = 10000,
                                                K::TU = FixedPoint(0.9),
                                                tol::FT = 1e-4,
                                                tosave::Vector{String} = String[],
                                                force_dict = Dict()) where {TM <: FPModel, TI <: IntegrationMethod,  TU<:Solver}

    nOP = n_order_params(Model)
    nOconjP = n_conj_order_params(Model)
    out = zeros(length(O0vec), 2*nOP  + 3 + 1 + nOconjP + length(tosave)) # that plus 1 is the ranking of free energies, from the most dominant to the least
    for (iO0,O0) in enumerate(O0vec)
        Model.O .= deepcopy(O0)
        conv, iter, vareps = solve_SPequations!(Model,X;niter = niter, K = K, tol = tol, force_dict = force_dict,showinfo= false,printprogress = false)
        obs_tosave = copy([Model.aux[oo] for oo in tosave])
        out[iO0,:] .= deepcopy([O0; Model.O; Model.aux["f"]; conv; iter; -12345.0;  get_Oconj(Model); obs_tosave])
    end

    # now fill the column with 1.0 with the ranking
    all_free_energies = out[:,2*nOP + 1]
    @assert unique(out[:, 2*nOP + 4]) == [-12345.0]

    out[:, 2*nOP + 4] .= denserank(round.(all_free_energies, digits =  log10(tol) +1 |> abs |> floor |> Int ))

    # check this 
    for (iO0,O0) in enumerate(O0vec)
        if all_free_energies[iO0] == NaN
            out[iO0, 2*nOP + 4] = length(O0vec)
        end
    end
    return out
end

