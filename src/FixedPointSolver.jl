using ExtractMacro, FastGaussQuadrature, LinearAlgebra, Statistics, ProgressMeter, StatsBase, Distributions, OrderedCollections 
using NamedArrays, Printf, ValueHistories, SpecialFunctions, Intervals, BigCombinatorics, MultinomialSeries
using DelimitedFiles, InteractiveUtils, Roots, OffsetArrays #UniformIsingModels
using QuadGK, Trapz, Cubature, HCubature, PrettyTables, Dates


include("types.jl")
include("numerical_integration.jl")
include("solvers.jl")
include("iterate.jl")
include("utils.jl")
include("greedy_minimum.jl")
include("scanning.jl")
include("criticalpoints.jl")
include("parse.jl")

include("discrete_expectation.jl")

# traces 
include("traces/HomogeneousSpin1half.jl")
include("traces/HomogeneousSpin1.jl")
include("traces/Spin1half.jl")

include("models/models_list.jl")

model_list = vcat(subtypes.(subtypes(FPModel))...) ∪ subtypes(FPModel) |> unique
;

critical_functions = collect_functions_with_underscore(Main);
; 



# export FPModel,DiscreteExpectation,IntegrationMethod
# for subtype in model_list
#     type_name = Base.typename(subtype)
#     @eval export $(Symbol(type_name.name))
# end
# for subtype in subtypes(IntegrationMethod)
#     type_name = Base.typename(subtype)
#     @eval export $(Symbol(type_name.name))
# end

# export create_model,lhs!,free_energy, compute_aux!
# export solve_SPequations!,solve_Multiple_initialconditions
# export Scan, ScannedModel, scanning
# export ScanCritical,ScannedCriticalModel,findcriticalpoint


