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
    "--OP"
        nargs = 2 # first argument name, second argument value, type inferred from the model type
        action = :append_arg
    "--parametric_type"
        nargs = '*'
        default = ["";""]
####### parameters to loop over #######
    "--param_loop"
        nargs = '*'
####### parameters scan over #######    
    "--pscan"
        arg_type = String
        required = true
        ####### SCANNING PARAMETERS
    "--pstart"
        arg_type = Float64
        default = 0.1
    "--delta_p"
        arg_type = Float64
        default = 0.1
    "--p_extrema"
        nargs = 2
        arg_type = Float64
        default = [-Inf; Inf]
    "--protocol"
        arg_type = String    
        default = "follow"
    "--tol_scancrit"
        arg_type = Float64
        default = 1e-4
    "--increase"
        arg_type = Bool
        default = true
    "--critical_line_name"
        arg_type = String
        default = "Prova"
    "--function"
        arg_type = String
#        required = true
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
    "--dampfirst"
        arg_type = Float64
        default = 0.96
####### MISC  #######
    "--aux_tosave"
        nargs = '*'
        arg_type = String
        default = String[]
    "--saveeachscan"
        arg_type = Bool    
        default = true
    "--dir"
        arg_type = String
        default = "res/prove"
end

p_args = parse_args(parse_settings)

dir = p_args["dir"]
create_subdir(dir,"Scans")


############# Get Model type, names and fields of its control parameters #############
ModelType, cparams_names, cparams_types = parse_Model_type(p_args["Model"])
nparams = length(cparams_names)

pscan = parse_inputsymbols(p_args["pscan"])
ploop = parse_inputsymbols(p_args["param_loop"][1])



#### fixed_params is the list of the ones you are giving to the model
fixed_params, vec_params = parse_model_parameters(p_args["param"],cparams_names, cparams_types, [pscan;ploop])


tuple_params = merge(vec_params...)
Model = set_startingModel(ModelType,tuple_params,p_args["parametric_type"])

@assert ploop in cparams_names
idx_ploop = findfirst(x->x==ploop,cparams_names)

ploop_type = cparams_types[findfirst(x->x==ploop,cparams_names)]
ploop_values = parse_inputarg_param_loop(p_args["param_loop"],ModelType)[2]

function_name = p_args["function"]

function singlerun(ModelType,ploop,ploop_v,pscan,tuple_params, p_args, critical_line_log_file)
    thid = Threads.threadid()
    
    X = parse_integration_method(p_args)
    
    
    Model = create_model(set_startingModel(ModelType,tuple_params,p_args["parametric_type"]), reset_parameter(tuple_params,ploop,ploop_v) )
    # set order parameters
    op_input = collect([parse_inputsymbols(x[1]) for x in p_args["OP"]])
    for n in names(Model.O,1)
        idx_op = findfirst(x->x==n,op_input)
        Model.O[n] = parse(Float64,p_args["OP"][idx_op][2])
    end
    
    Ostart = deepcopy(Model.O)
    scancrit = ScanCritical(pscan,p_args["pstart"], p_args["delta_p"], p_args["tol_scancrit"], 0.5,p_args["p_extrema"],Symbol(p_args["protocol"]),p_args["increase"])
    upordown = p_args["increase"] ? "increase" : "decrease"
    filepath = string(p_args["critical_line_name"],"_", string(Model),"_",replace_Greek(params_to_textstring(Model;exclude = [pscan])),"_scan_over_",replace_Greek(string(pscan)),"_",upordown,"_loop_over_",replace_Greek(string(ploop)),print_integration_params(X),"_protocol_",scancrit.protocol, "_tol_",p_args["tol"],"_tolc_",scancrit.Δptol)
    
    teps = @elapsed out = findcriticalpoint(Model,scancrit, critical_functions[p_args["function"]],X;
                                        tol = p_args["tol"], K=FixedPoint(p_args["damp"]),Kfirstrun=FixedPoint(p_args["dampfirst"]), niter = p_args["maxiter"], 
                                        printprogress = false,fsave = string(dir,"/Scans/",filepath, ".txt"),
                                        tosave = p_args["aux_tosave"]  )

    
    writedlm(string(dir,"/Scans/",filepath, "_RES.txt"), [ploop_v; out.pfinal; out.pfinallastconv; teps] )

    log_message(string("on thread ",thid, "Line -> ",p_args["critical_line_name"]," on Model = ",string(ModelType)," -> params  ",print_tuple_params(Model; exclude = [pscan]), " --> critical ", pscan,"= ", (out.pfinal, out.pfinallastconv), " Init cond.  ",Ostart, " elapsed=",teps," seconds" ),critical_line_log_file)
end

starting_info = """
******************************* CRITICAL LINE ******************************* 
Model: $(string(Model; separator = " "))
Fixed parameters: $(params_to_textstring(tuple_params;separator="; ",exclude=[ploop;pscan],equalsymbol="="))
Critical line name $(p_args["critical_line_name"])
Loop parameters: $ploop,  number of parallel scan is $(length(ploop_values)) --- Extrema = $(extrema(ploop_values))
scan parameter: $pscan ###  protocol = $(p_args["protocol"])   --- function $function_name  
*****************************************************************************
"""
print_and_log(starting_info,critical_line_log_file)


################# set up loop ##################
number_jobs = length(ploop_values)
job_idx = collect(1:number_jobs);

const completed_jobs = Threads.Atomic{Int}(0)

##### 
using Random;  Random.shuffle!(job_idx) # permute values so to avoid that slower-convergence runs to be assigned to same thread
#####

progress=set_progress_bar(number_jobs,"critical_line")
timestart = Dates.now()

Threads.@threads for iploop in job_idx
    singlerun(ModelType,ploop,ploop_values[iploop],pscan,tuple_params,p_args,critical_line_log_file)
    # Update the progress bar with a message showing completed jobs
    Threads.atomic_add!(completed_jobs, 1)
    next!(progress; showvalues = [(:jobs_completed, "$(completed_jobs[]) / $number_jobs" )] )
end

timeend = Dates.now()

elapsed = canonicalize(round(timeend-timestart, Second(1)))
print_and_log(string("TOTAL TIME ELAPSED: ", elapsed), critical_line_log_file)

########################################################################################################################
####################################################### WRAP DATA ######################################################
########################################################################################################################

values_final = zeros(length(ploop_values),3)
teps_stats = zeros(length(ploop_values))
scans = Dict{ploop_type,Matrix}()


X = parse_integration_method(p_args)
upordown = p_args["increase"] ? "increase" : "decrease"

single_file = zeros(4)
progress = Progress(length(ploop_values), desc="Wrapping", barlen=settings["barlen"],color = progress_color["critical_line"], barglyphs=BarGlyphs("[=> ]"))
for iploop in eachindex(ploop_values)
    ploop_v = ploop_values[iploop]
    params = (; tuple_params..., ploop=>ploop_v)
    filepath = string(p_args["critical_line_name"],"_", string(Model),"_",replace_Greek(params_to_textstring(params;exclude = [pscan])),"_scan_over_",replace_Greek(string(pscan)),"_",upordown,"_loop_over_",replace_Greek(string(ploop)),print_integration_params(X),"_protocol_",p_args["protocol"], "_tol_",p_args["tol"],"_tolc_",p_args["tol_scancrit"])

    single_file .= readdlm(string(dir,"/Scans/",filepath, "_RES.txt"))
    values_final[iploop,:] .=  deepcopy(single_file[1:3])
    teps_stats[iploop] = single_file[4]
    rm(string(dir,"/Scans/",filepath, "_RES.txt"))

    if p_args["saveeachscan"]
        scans[ploop_v] = readdlm(string(dir,"/Scans/",filepath, ".txt"))
    end
    rm(string(dir,"/Scans/",filepath, ".txt"))
    next!(progress)
end

final_path = string(p_args["critical_line_name"],"_", string(Model),"_",replace_Greek(params_to_textstring(tuple_params;exclude = [ploop; pscan])),"_scan_over_",replace_Greek(string(pscan)),"_",upordown,"_loop_over_",replace_Greek(string(ploop)),print_integration_params(X),"_protocol_",p_args["protocol"], "_tol_",p_args["tol"],"_tolc_",p_args["tol_scancrit"])
print_and_log(string("filename: ",final_path), critical_line_log_file)
print_and_log(string("directory: ",dir), critical_line_log_file)

save(string(dir, "/",final_path,".jld"), "input_args", p_args, 
                                        "values_final", values_final,
                                        "ploop_values",ploop_values,
                                        "scans", scans, 
                                        "filename",final_path,
                                        "time_start",timestart, 
                                        "time_end", timeend,
                                        "teps_stats",teps_stats, 
                                        "elapsed",elapsed,
                                        "host",get_hostname(),                                        
                                        "git_info",get_commit_info()
)

end_log(critical_line_log_file)