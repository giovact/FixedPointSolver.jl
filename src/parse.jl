function parse_Model_type(str_model::String)
    ModelType = model_list[findfirst(string(TM) == p_args["Model"] for TM in model_list)]
    cparams_names = fieldtype(ModelType,:params).parameters[1]
    cparams_types = fieldtypes(fieldtype(ModelType,:params).parameters[2])
    return ModelType, cparams_names, cparams_types
end

function parse_model_parameters(param_list_toparse, cparams_names, cparams_types,set_dumb)
    fixed_params = collect([parse_inputsymbols(param_list_toparse[ip][1]) for ip=1:length(param_list_toparse)] )
    vec_params = []
    for (ip,p) in enumerate(cparams_names)
        if p ∈ set_dumb
            push!(vec_params,(;p=>dumb(cparams_types[ip])))
        else
            iip = findfirst(x->x==p,fixed_params)
            push!(vec_params,(;p=>parse(cparams_types[ip], param_list_toparse[iip][2] )))
        end
    end
    return fixed_params, vec_params
end


function parse_inputsymbols(s::String)
    has_a_greek = findall(isequal('%'), s)
    if isempty(has_a_greek)
        # means there is no greek letter
        return Symbol(s)
    else
        @assert length(has_a_greek)==2 # for now, might be that you have two greek letters in the same symbol but ok
        #println(s[1:has_a_greek[1]-1])
        #println(Symbol(latin_to_greek[s[has_a_greek[1]+1:has_a_greek[2]-1]]))
        #println(Symbol(s[has_a_greek[2]+1:end]))
        return Symbol(s[1:has_a_greek[1]-1],latin_to_greek[s[has_a_greek[1]+1:has_a_greek[2]-1]],s[has_a_greek[2]+1:end])
    end
end


function set_startingModel(modeltype, tuple_params, parametric)
    if isempty(parametric[1]) && isempty(parametric[2])
        return modeltype(tuple_params...)
    elseif parametric[1] == "Potential"
        subtypes_PotentialPerceptron = subtypes(Potential)
        LRtype = subtypes_PotentialPerceptron[findfirst(string(pot) == parametric[2] for pot in subtypes_PotentialPerceptron)]
        learning_rule = LRtype(parse(Bool,parametric[3]))
        return modeltype(tuple_params...,learning_rule)
    elseif parametric[1]=="DiscreteExpectation"
        if parametric[2]=="PureStatesRetrievalBinary"
            return modeltype(tuple_params...,PureStatesRetrievalBinary())
        elseif parametric[2]=="MixedStatesRetrievalBinary"
            return modeltype(tuple_params...,MixedStatesRetrievalBinary(parse(Int,parametric[3])))
        elseif parametric[2]=="MixedStatesRetrievalBinarySYM"
            return modeltype(tuple_params...,MixedStatesRetrievalBinarySYM(parse(Int,parametric[3])))
        end

        if parametric[2]=="PureStatesRetrievalSparsePatterns"
            return modeltype(tuple_params...,PureStatesRetrievalSparsePatterns(parse(FT,parametric[3])))
        elseif parametric[2]=="MixedStatesRetrievalSparsePatternsBinary"
            return modeltype(tuple_params...,MixedStatesRetrievalSparsePatternsBinary(parse(FT,parametric[3]), parse(Int,parametric[4])))
        elseif parametric[2]=="MixedStatesRetrievalSparsePatternsBinarySYM"
            return modeltype(tuple_params...,MixedStatesRetrievalSparsePatternsBinarySYM(parse(FT,parametric[3]), parse(Int,parametric[4])))
        end
        throw("Not supported yet")
    else
        throw("Not supported yet")
    end 
end


function set_startingModel(modeltype, parametric)
    if isempty(parametric[1]) && isempty(parametric[2])
        return modeltype()
    elseif parametric[1] == "Potential"
        subtypes_PotentialPerceptron = subtypes(Potential)
        LRtype = subtypes_PotentialPerceptron[findfirst(string(pot) == parametric[2] for pot in subtypes_PotentialPerceptron)]
        learning_rule = LRtype(parse(Bool,parametric[3]))
        return modeltype(learning_rule)
    elseif parametric[1]=="DiscreteExpectation"
        if parametric[2]=="PureStatesRetrievalBinary"
            modeltype(PureStatesRetrievalBinary())
        elseif parametric[2]=="MixedStatesRetrievalBinary"
            modeltype(MixedStatesRetrievalBinary(parse(Int,parametric[3])))
        elseif parametric[2]=="MixedStatesRetrievalBinarySYM"
            modeltype(MixedStatesRetrievalBinarySYM(parse(Int,parametric[3])))
        else
            throw("Not supported yet")
        end
    else
        throw("Not supported yet")
    end 
end


function parse_inputarg_param_loop_old(input, Modeltype)
    # the second output argument here is not necessarily sorted
    @assert length(input) == 5
    ploop = parse_inputsymbols(input[1])

    cparams_names = fieldtype(ModelType,:params).parameters[1]
    cparams_types = fieldtypes(fieldtype(ModelType,:params).parameters[2])
    plt = cparams_types[findfirst(x->x==ploop,cparams_names)]

    pstart = parse(plt, input[2])
    pstep = parse(plt, input[3])
    pstop = parse(plt, input[4])    
    
    
    if input[5] == "linear"
        return ploop, collect(pstart:pstep:pstop)  
    elseif input[5] == "reciprocal"
        return  ploop, sort(1 ./ collect(pstart:pstep:pstop))
    elseif input[5] == "exp10"
        return ploop, 10 .^ collect(pstart:pstep:pstop)
    else 
        throw("Choose a correct scaling of the parameters -> available now [linear; reciprocal; exp10]")
    end
    # instead of returning just the value, you couldreturn also the label (needed for visualization) and the corresponding scale you want for the plot (i.e. linear, or log)
    
end


function parse_inputarg_param_loop(input, Modeltype)
    # the second output argument here is not necessarily sorted
    #@assert length(input) == 6 || length(input) == 8
    ploop = parse_inputsymbols(input[1])

    cparams_names = fieldtype(ModelType,:params).parameters[1]
    cparams_types = fieldtypes(fieldtype(ModelType,:params).parameters[2])
    plt = cparams_types[findfirst(x->x==ploop,cparams_names)]

    ### handle list of specific values 
    if input[2] == "list"
        values = zeros(plt,0)
        for i = 3 : length(input)
            push!(values, parse(plt, input[i]))
        end 
        return ploop, values
    end

    pstart = parse(plt, input[2])
    pstop = parse(plt, input[3])   

    @assert pstop >= pstart

    pstep = dumb(plt)
    mesh_length=0
    if input[4]=="step"
        pstep=parse(plt,input[5])
    elseif input[4]=="length"
        mesh_length=parse(Int,input[5])
    else
        throw("invivalid option for mesh -> choose among [step, ploop_type] or [length, Int]")
    end

    ## define mesh 
    mesh = mesh_length>0 ? collect(range(pstart,stop = pstop, length=mesh_length)) : collect(range(pstart,stop = pstop, step=pstep))
    if length(input) == 8
         if input[7] == "scaleup"
            mesh .*= parse(plt,input[8])
         elseif input[7] == "scaledown"
            mesh ./= parse(plt,input[8])
         else
            throw("select valid scale option -> choose among [scaleup, scaledown]") 
        end
    end

    if input[6] == "linear"
        return ploop, mesh
    elseif input[6] == "reciprocal"
        return  ploop, sort(1 ./ mesh )
    elseif input[6] == "exp10"
        return ploop, 10 .^ mesh
    else 
        throw("Choose a correct scaling of the parameters -> available now [linear; reciprocal; exp10]")
    end
    # instead of returning just the value, you couldreturn also the label (needed for visualization) and the corresponding scale you want for the plot (i.e. linear, or log)
    
end


function get_hostname()
    if "HOSTNAME" in keys(ENV)
        return ENV["HOSTNAME"]
    elseif "SESSION_MANAGER" in keys(ENV)
        return ENV["SESSION_MANAGER"]
    else
        return "Not available"
    end
end


parse_integration_method(args) = set_integration_method(args["IntMethod"], args["npointsInt"], args["bound"])

# solver from the command line (--solver, --anderson_m); damp is ρ for both solvers (same convention, see solvers.jl)
function parse_solver(args, damp::Float64)
    args["solver"] == "FixedPoint" && return FixedPoint(damp)
    args["solver"] == "Anderson"   && return Anderson(args["anderson_m"]; ρ = damp)
    throw("Choose a valid solver -> available now [FixedPoint; Anderson]")
end