struct Scan
    param::Symbol
    pvalues::Vector
    protocol::Symbol
    reverseordering::Bool
    function Scan(param, pvalues,protocol, reverseordering)
        @assert issorted(pvalues)
        protocol == :restart && reverseordering && @warn "Useless to use reverse ordering when restarting"
        protocol ∉ [:follow; :restart] && throw("select valid protocol") 
        new(param, pvalues,protocol, reverseordering)
    end
end

Scan(p,v,protocol; reverse = false) = Scan(p,v,protocol, reverse)

struct ScannedModel{TM<:FPModel}
    Model::TM
    scan::Scan
    data::Matrix{FT}
    names::Vector{String}
end


function scanning_info(scan::Scan, Modelstart)
    scan.protocol == :restart && return string(" RESTART - each time from -> ",print_op(Modelstart) )
    scan.protocol == :follow && return string(" FOLLOW fps; reverse=", scan.reverseordering)
end

function scanning(Modelstart::TM,scan::Scan,X::TI;
                                niter::Int = 10000,
                                K::TU = FixedPoint(0.9),
                                Kfirstrun::TU = FixedPoint(0.95),
                                tol::FT = 1e-5,
                                eps_reset::FT = -Inf,
                                printprogress::Bool = true,
                                tosave::Vector{String} = String[],
                                fsave::String = "",
                                force_dict = Dict{Symbol, FT}(),
                                overridetype::Bool = false,
                                ) where {TM<:FPModel, TI<:IntegrationMethod, TU<:Solver }
    
    @extract scan : param pvalues protocol reverseordering

    pvalues_sorted  =  reverseordering ? reverse(pvalues) : pvalues
    
    O0 = deepcopy(Modelstart.O)
    Model = create_model(deepcopy(Modelstart), reset_parameter(Modelstart.params, param,overridetype ? Modelstart.params[param] : -Inf ) ) # dummy. Such that Model.O  == Model.Ostart
    @assert typeof(Model) == typeof(Modelstart) #debug
    @assert Model.O == Modelstart.O

    infotoprint = string("Scanning over parameter ",param, "  ", print_params(Model;separator= " - ",exclude =[param]), scanning_info(scan, Modelstart))
    prog = printprogress ? Progress( length(pvalues), infotoprint ) : nothing

    data = zeros( length(pvalues), n_order_params(Modelstart) + 1 + 2 + n_conj_order_params(Modelstart) + length(tosave) + 2 ) # O, f, conv, iter, Oconj if present, then aux to save, then elapsed time in seconds
    
    names_keys = [string.(names(Modelstart.O,1)); "f"; "conv";"iter";]
    if has_conj_order_params(Modelstart)
        names_keys = [names_keys; get_Oconjnames(Modelstart)]
    end
    if !isempty(tosave)
        names_keys = [names_keys;tosave]
    end
    names_keys = [names_keys; "elapsed";"elapsed_cumulated"]
    teps_cum = 0.0    
    for (i, v) in enumerate(pvalues_sorted)

        Model = create_model(Model, reset_parameter(Model.params, param,v) )
        teps = @elapsed conv, iter, vareps = solve_SPequations!(Model,X; K = scan.protocol == :restart ? K : Kfirstrun, tol = tol, niter = niter,showinfo=false,printprogress=false,force_dict = force_dict)
        obs_tosave = copy([Model.aux[oo] for oo in tosave])
        
        teps_cum +=teps
        data[i,:] .= copy([Vector(Model.O); Model.aux["f"]; conv; iter; get_Oconj(Model); obs_tosave; teps; teps_cum]) # probably ValueHistories is cooler here
        
        (length(fsave)>0) && writedlm(fsave,hcat(pvalues_sorted, data))
        
        if protocol == :restart
            Model.O .= deepcopy(O0)
        elseif protocol == :follow 
            reset_orderparams!(Model.O,eps_reset)
        end

        #printprogress && ProgressMeter.next!(prog; showvalues = [print_op_withparam(Model,param); (:f,Model.aux["f"]); (:iter,iter); (:conv, conv)] )
        printprogress && ProgressMeter.next!(prog; showvalues = print_scanning(Model, param, iter, conv, vareps, tosave) )

    end
    
    if reverseordering
        return ScannedModel(Modelstart,scan,reverse(data, dims = 1), names_keys)
    else
        return ScannedModel(Modelstart,scan, data, names_keys)
    end
end


function print_scanning(Model::TM,param::Symbol, iter::Int, conv::Int,vareps::FT, tosave::Vector{String}) where {TM <: FPModel}
    
    isempty(tosave) && return [print_op_withparam(Model,param); (:f,Model.aux["f"]); (:iter,iter); (:conv, conv); (:ε, vareps)]
    vcat([print_op_withparam(Model,param); (:f,Model.aux["f"]); (:iter,iter); (:conv, conv); (:ε, vareps); ], [(Symbol(x),Model.aux[x]) for x in tosave])
end


## highlight

function showtable(S::ScannedModel) 
    data_complete = hcat(S.scan.pvalues, S.data)
    hl_1 = Highlighter(f = (data_complete,i,j) -> data_complete[i,3+n_order_params(S.Model)]==1, crayon = crayon"green")
    hl_0 = Highlighter(f = (data_complete,i,j) -> data_complete[i,3+n_order_params(S.Model)]==0, crayon = crayon"yellow")
    hl_m1 = Highlighter(f = (data_complete,i,j) -> data_complete[i,3+n_order_params(S.Model)]==-1, crayon = crayon"red")

    pretty_table(data_complete, header = [string(S.scan.param); S.names],  header_crayon = crayon"blue bold", highlighters = (hl_1,hl_0,hl_m1))
end

