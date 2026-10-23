include("../src/FixedPointSolver.jl")
include("settings.jl")
include("logging.jl")

using DelimitedFiles, ArgParse, JLD, InteractiveUtils

parse_settings = ArgParseSettings()

@add_arg_table! parse_settings begin
####### Model and fixed parameters #######
    "--Model"
        arg_type = String
        required = true
    "--param"
        nargs = 2 # first argument name, second argument value, type inferred from the model type
        action = :append_arg
    "--parametric_type"
        nargs = '*'
        default = ["";"" ]        
####### parameters to loop over #######
    "--param_loop"
        nargs = '*'
        required = true
        action = :append_arg
    "--OPnames"
        arg_type = String
        nargs = '*' # just names
        required = true
    "--OP0"
        nargs = '*' # first argument name, second argument value, type inferred from the model type
        arg_type = Float64
        action = :append_arg
####### parameters for solving each SP equation
    "--IntMethod"
        arg_type = String
        default = "LegendreQuadrature"
    "--npointsInt"
        arg_type = Int
        default = 200
    "--bound"
        arg_type = Float64
        default = 6.0
    "--tol"
        arg_type = Float64
        default = 1e-6
    "--maxiter"
        arg_type = Int
        default = 3000
    "--damp"
        arg_type = Float64
        default = 0.92
    "--solver"
        arg_type = String
        default = "FixedPoint"
        range_tester = x -> x in ["FixedPoint", "Anderson"]
    "--anderson_m"
        arg_type = Int
        default = 5
####### MISC  #######
    "--suffix" # if needed to distinguish two files created with the same settings 
        arg_type = String
        default = ""
    "--aux_tosave"
        nargs = '*'
        arg_type = String
        default = String[]
    "--dir"
        arg_type = String
        default = "res/prove"
    "--plot"
        arg_type=Bool
        default=true
    "--overwrite"
        arg_type=Bool
        default=false
end

p_args = parse_args(parse_settings)

dir = p_args["dir"]
create_subdir(dir,"SingleRuns")


############# Get Model type, names and fields of its control parameters #############
ModelType, cparams_names, cparams_types = parse_Model_type(p_args["Model"])
nparams = length(cparams_names)


loop_params = collect([parse_inputsymbols(p_args["param_loop"][ip][1]) for ip=1:length(p_args["param_loop"])] )
@assert length(loop_params)==2
fixed_params, vec_params = parse_model_parameters(p_args["param"],cparams_names, cparams_types, loop_params)

ploop_values = Dict(loop_params[ip]=>parse_inputarg_param_loop(p_args["param_loop"][ip], ModelType)[2] for ip in eachindex(p_args["param_loop"]) )

@assert length(fixed_params) + length(loop_params) == nparams

tuple_params = merge(vec_params...)
Model = set_startingModel(ModelType,tuple_params,p_args["parametric_type"])

starting_info = """
*************** EXPLORATION with different initial conditions *************** 
Model: $(string(Model; separator = " "))
Fixed parameters: $(params_to_textstring(tuple_params;separator="; ",exclude=loop_params,equalsymbol="="))
Loop parameters: $loop_params
Total number of mesh points in $(loop_params[1]) = $(length(ploop_values[loop_params[1]]))  Extrema = $(extrema(ploop_values[loop_params[1]]))
Total number of mesh points in $(loop_params[2]) = $(length(ploop_values[loop_params[2]]))  Extrema = $(extrema(ploop_values[loop_params[2]]))
Total number of runs = $(length(ploop_values[loop_params[1]]) * length(ploop_values[loop_params[2]]))
Solver: $(parse_solver(p_args, p_args["damp"]))
*****************************************************************************
"""
print_and_log(starting_info,exploration_log_file)



function singlerun(ModelType,loop_params,pl_values,tuple_params, p_args,exploration_log_file)
    thid = Threads.threadid()
    
    X = parse_integration_method(p_args)
    
    Model = create_model(set_startingModel(ModelType,tuple_params,p_args["parametric_type"]), reset_parameter(reset_parameter(tuple_params,loop_params[1],pl_values[1]),loop_params[2],pl_values[2]) )

    # set order parameters, for the moment, assert equal ordering between names and Model
    @assert length(p_args["OPnames"]) == n_order_params(Model)
    for (in,n) in enumerate(names(Model.O,1))
        @assert n==parse_inputsymbols(p_args["OPnames"][in])
    end
    
    filepath = string("MultipleInit_", string(Model),"_",suffix(p_args["suffix"]),replace_Greek(params_to_textstring(Model)),print_integration_params(X),"_tol_",p_args["tol"])
    teps = @elapsed out = solve_Multiple_initialconditions(Model,X,p_args["OP0"],tol = p_args["tol"], K=parse_solver(p_args, p_args["damp"]), niter = p_args["maxiter"], tosave = p_args["aux_tosave"]  )
    writedlm(string(dir,"/SingleRuns/",filepath, ".txt"), out )

    log_message(string("on thread ",thid, " Model -> ",ModelType,"  --- params=  ",print_tuple_params(Model), " --> done  in ", teps,"  seconds"), exploration_log_file)
end


#### check if file exists 
X = parse_integration_method(p_args)
final_path = string("MultipleInit_", string(Model),"_",suffix(p_args["suffix"]),replace_Greek(params_to_textstring(tuple_params;exclude = loop_params)),print_integration_params(X),"_tol_",p_args["tol"])

if isfile(string(dir, "/",final_path,".jld")) && p_args["overwrite"]==false
    print_and_log("FILE ALREADY PRESENT, skipping",exploration_log_file)
    print_and_log(string("filename  --> ",final_path), exploration_log_file)
    print_and_log(string("directory --> ",dir), exploration_log_file)
    if p_args["plot"]
        print_and_log("NOW PLOTTING", exploration_log_file)
        include("../src/visualization/Plotheaders.jl")
        plot_res_exploration(string(dir, "/",final_path))
    end
    exit()
end    


# indices -> 
Cidx = CartesianIndices(( length(ploop_values[loop_params[1]]),length(ploop_values[loop_params[2]])) ) |> collect |> vec
tuples_loop = [(ploop_values[loop_params[1]][CI.I[1]],ploop_values[loop_params[2]][CI.I[2]]) for CI in Cidx]

################# set up loop ##################
number_jobs = length(Cidx)
job_idx = collect(1:number_jobs);

const completed_jobs = Threads.Atomic{Int}(0)

##### 
using Random;  Random.shuffle!(job_idx) # permute values so to avoid that slower-convergence runs to be assigned to same thread
#####

progress=set_progress_bar(number_jobs,"exploration")
timestart = Dates.now()

Threads.@threads for idx in job_idx
    singlerun(ModelType,loop_params,tuples_loop[idx],tuple_params,p_args,exploration_log_file)
    # Update the progress bar with a message showing completed jobs
    Threads.atomic_add!(completed_jobs, 1)
    next!(progress; showvalues = [(:jobs_completed, "$(completed_jobs[]) / $number_jobs" )] )
end

timeend = Dates.now()
elapsed = canonicalize(round(timeend-timestart, Second(1)))

print_and_log(string("TOTAL TIME ELAPSED: ", elapsed), exploration_log_file)


########################################################################################################################
####################################################### WRAP DATA ######################################################
########################################################################################################################


# output in the form [O0; Model.O; Model.aux["f"]; conv; iter; get_Oconj(Model); obs_tosave]
# la cosa più semplice forse è salvarti un array a 3 dimensioni per ogni valore salvato, dizionario di array a 3 dimensioni
out_keys = vcat(String.(names(Model.O,1)), ["f";"conv";"iter";"f_ranking"])
n_conj_order_params(Model)!=0 && append!(out_keys,get_Oconjnames(Model) )
!isempty(p_args["aux_tosave"])!=0 && append!(out_keys,p_args["aux_tosave"])


wrapped_out = Dict(k=>zeros(length(ploop_values[loop_params[1]]),length(ploop_values[loop_params[2]]),length(p_args["OP0"])) for k in out_keys)
nop = n_order_params(Model) # index to start from

progress = Progress(length(Cidx), desc="Wrapping", barlen=settings["barlen"],color = progress_color["exploration"], barglyphs=BarGlyphs("[=> ]"))
for (idx, CI) in enumerate(Cidx)

    
    params = (; tuple_params..., loop_params[1]=>tuples_loop[idx][1],loop_params[2]=>tuples_loop[idx][2] )

    #check on CartesianIndices
    @assert params[loop_params[1]] == ploop_values[loop_params[1]][CI.I[1]]
    @assert params[loop_params[2]] == ploop_values[loop_params[2]][CI.I[2]] 
    filepath = string("MultipleInit_", string(Model),"_",suffix(p_args["suffix"]),replace_Greek(params_to_textstring(params)),print_integration_params(X),"_tol_",p_args["tol"])
    res = readdlm(string(dir,"/SingleRuns/",filepath, ".txt") )
    @assert n_order_params(Model) + length(out_keys) == size(res,2)
    rm(string(dir,"/SingleRuns/",filepath, ".txt"))
    # now wrap the dictionary


    for (in,n) in enumerate(names(Model.O,1))
        wrapped_out[String(n)][CI.I[1], CI.I[2], :] .= res[:,nop+in]
    end
    wrapped_out["f"][CI.I[1], CI.I[2], :] .= res[:,2*nop+1]
    wrapped_out["conv"][CI.I[1], CI.I[2], :] .= res[:,2*nop+2]    
    wrapped_out["iter"][CI.I[1], CI.I[2], :] .= res[:,2*nop+3]
    wrapped_out["f_ranking"][CI.I[1], CI.I[2], :] .= res[:,2*nop+4]

    if has_conj_order_params(Model)
        for (in,n) in enumerate(names(Model.Oconj,1))
            wrapped_out[string(n, "conj")][CI.I[1], CI.I[2], :] .= res[:,2 * nop + 4 + in]
        end
    end
    if !isempty(p_args["aux_tosave"])
        for (in,n) in enumerate(p_args["aux_tosave"])
            wrapped_out[n][CI.I[1], CI.I[2], :] .= res[:,2 * nop + 4 + n_conj_order_params(Model)+in]
        end
    end
    next!(progress)
end

print_and_log(string("filename: ",final_path), exploration_log_file)
print_and_log(string("directory: ",dir), exploration_log_file)

save(string(dir, "/",final_path,".jld"), "input_args", p_args, 
                                        "loop_params", loop_params, 
                                        "ploop_values", ploop_values, 
                                        "OPnames",String.(names(Model.O,1)),
                                        "res",wrapped_out , 
                                        "filename",final_path,
                                        "time_start",timestart, 
                                        "time_end", timeend, 
                                        "elapsed",elapsed,
                                        "host",get_hostname(),
                                        "git_info",get_commit_info()
)


end_log(exploration_log_file) 

if p_args["plot"]
    print_and_log("NOW PLOTTING", exploration_log_file)
    include("../src/visualization/Plotheaders.jl")
    plot_res_exploration(string(dir, "/",final_path))
end



