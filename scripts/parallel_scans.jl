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
####### parameter to loop over #######
    "--param_loop"
        nargs = '*'
####### parameters scan over #######    
    "--param_scan"
        nargs = '*'
        required = true
    "--protocol"
        arg_type = String    
        default = "follow"
    "--increase"
        arg_type = Bool
        default = true
##### initial condition
    "--OPnames"
        arg_type = String
        nargs = '*' # just names
        required = true
    "--OP0"
        nargs = '*'
        arg_type = Float64
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
        
end

p_args = parse_args(parse_settings)

dir = p_args["dir"]
create_subdir(dir,"Scans")

############# Get Model type, names and fields of its control parameters #############
ModelType, cparams_names, cparams_types = parse_Model_type(p_args["Model"])
nparams = length(cparams_names)

ploop, ploop_values = parse_inputarg_param_loop(p_args["param_loop"],ModelType)
pscan, pscan_values = parse_inputarg_param_loop(p_args["param_scan"],ModelType)

fixed_params, vec_params = parse_model_parameters(p_args["param"],cparams_names, cparams_types, [ploop; pscan])

tuple_params = merge(vec_params...)
Model = set_startingModel(ModelType,tuple_params,p_args["parametric_type"])

scan = Scan(pscan,sort(pscan_values), Symbol(p_args["protocol"]), reverse = !p_args["increase"])

function singlerun(ModelType,ploop,ploop_v,scan, tuple_params,p_args,parallel_scans_log_file)
    thid = Threads.threadid()
    
    X = parse_integration_method(p_args)

    Model = create_model(set_startingModel(ModelType,tuple_params,p_args["parametric_type"]), reset_parameter(tuple_params,ploop,ploop_v) )
    # set order parameters

    @assert length(p_args["OPnames"]) == n_order_params(Model)
    for (in,n) in enumerate(names(Model.O,1))
        @assert n==parse_inputsymbols(p_args["OPnames"][in])
        Model.O[n] = p_args["OP0"][in]
    end
    
    Ostart = deepcopy(Model.O)
    
    filepath = string("Scanning_", string(Model),"_",suffix(p_args["suffix"]),replace_Greek(params_to_textstring(Model;exclude = [scan.param])),print_integration_params(X),"_tol_",p_args["tol"],"_protocol_",p_args["protocol"])

    teps = @elapsed out = scanning(Model,scan, X;
                                        tol = p_args["tol"], K=parse_solver(p_args, p_args["damp"]),Kfirstrun=parse_solver(p_args, p_args["dampfirst"]), niter = p_args["maxiter"], 
                                        printprogress = false,
                                        tosave = p_args["aux_tosave"],
                                        fsave = string(dir,"/Scans/",filepath, "_scan.txt")  )

    
    writedlm(string(dir,"/Scans/",filepath, "_scan.txt"), out.data)
    log_message(string("on thread ",thid, " Scanning -> over ", scan.param, " - extrema =", extrema(scan.pvalues), " reverse_ordering=", scan.reverseordering, "  on Model = ",string(ModelType)," -> params  ",print_tuple_params(Model; exclude = [scan.param]), "  init --> ",print_op(Ostart), " elapsed=",teps," seconds" ),parallel_scans_log_file)
    
end

starting_info = """
******************************* Parallel Scans ****************************** 
 Model: $(string(Model; separator = " "))
 fixed parameters: $(params_to_textstring(tuple_params;separator="; ",exclude=[ploop;pscan],equalsymbol="="))
 loop parameter: $ploop, number of parallel scan is $(length(ploop_values))
 scan parameter: $pscan, number of mesh points for each scan $(length(pscan_values)), protocol = $(p_args["protocol"])
 Solver: $(parse_solver(p_args, p_args["damp"])) (first run: $(parse_solver(p_args, p_args["dampfirst"])))
*****************************************************************************
"""
print_and_log(starting_info,parallel_scans_log_file)


################# set up loop ##################
number_jobs = length(ploop_values)
job_idx = collect(1:number_jobs);

const completed_jobs = Threads.Atomic{Int}(0)

##### 
using Random;  Random.shuffle!(job_idx) # permute values so to avoid that slower-convergence runs to be assigned to same thread
#####

progress=set_progress_bar(number_jobs,"parallel_scans")
timestart = Dates.now()

Threads.@threads for iploop in job_idx
    singlerun(ModelType,ploop,ploop_values[iploop],scan,tuple_params,p_args,parallel_scans_log_file)
   
    # Update the progress bar with a message showing completed jobs
    Threads.atomic_add!(completed_jobs, 1)
    next!(progress; showvalues = [(:jobs_completed, "$(completed_jobs[]) / $number_jobs" )] )end

timeend = Dates.now()
elapsed = canonicalize(round(timeend-timestart, Second(1)))

print_and_log(string("TOTAL TIME ELAPSED: ", elapsed), parallel_scans_log_file)

########################################################################################################################
####################################################### WRAP DATA ######################################################
########################################################################################################################


# dictionary for each value in ploop, add function to convert dictionary into matrix
scans_dict = Dict{typeof(ploop_values[1]),Matrix}()

X = parse_integration_method(p_args)

# do a dummy scan just to get names
dummy_scan_out = scanning(Model, scan,X; tosave = p_args["aux_tosave"], niter=1,printprogress=false);

progress = Progress(length(ploop_values), desc="Wrapping", barlen=settings["barlen"],color = progress_color["parallel_scans"], barglyphs=BarGlyphs("[=> ]"))
for iploop in eachindex(ploop_values)
    ploop_v = ploop_values[iploop]
    params = (; tuple_params..., ploop=>ploop_v)
    filepath = string("Scanning_", string(Model),"_",suffix(p_args["suffix"]),replace_Greek(params_to_textstring(params;exclude = [pscan])),print_integration_params(X),"_tol_",p_args["tol"],"_protocol_",p_args["protocol"])

    single_data=readdlm(string(dir,"/Scans/",filepath, "_scan.txt"))
    scans_dict[ploop_v] = single_data
    rm(string(dir,"/Scans/",filepath, "_scan.txt"))
    next!(progress)
end

final_path = string("ParallelScans_",string(Model),"_",suffix(p_args["suffix"]),replace_Greek(params_to_textstring(tuple_params;exclude = [ploop; pscan])),"_scan_over_",replace_Greek(string(pscan)),"_reverse_",!p_args["increase"],"_loop_over_",replace_Greek(string(ploop)),print_integration_params(X),"_protocol_",p_args["protocol"], "_tol_",p_args["tol"])
print_and_log(string("filename: ",final_path), parallel_scans_log_file)
print_and_log(string("directory: ",dir), parallel_scans_log_file)
save(string(dir, "/",final_path,".jld"), "input_args", p_args, 
                                        "pscan",pscan,
                                        "pscan_values",pscan_values,
                                        "ploop",ploop,
                                        "ploop_values",ploop_values,
                                        "scans", scans_dict,
                                        "names",dummy_scan_out.names,
                                        "filename",final_path,
                                        "time_start",timestart, 
                                        "time_end", timeend, 
                                        "elapsed",elapsed,
                                        "host",get_hostname(),
                                        "git_info",get_commit_info()
)


end_log(parallel_scans_log_file)
