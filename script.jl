import Pkg;
# Pkg.build("SpecialFunctions")
# Pkg.add("Gurobi")
# Pkg.add("JuMP")
# Pkg.add("RDatasets")
# Pkg.add("Plots")
# Pkg.add("DataTables")
# Pkg.add("Clustering")
# Pkg.add("LinearAlgebra")
# Pkg.add("Distributions")
# Pkg.add("StatsBase")
# Pkg.add("DelimitedFiles")
#Pkg.add("ThreadPools")
using RDatasets, DataTables, Clustering, LinearAlgebra, Distributions, StatsBase, DelimitedFiles, ThreadPools, Random
ENV["JULIA_NUM_THREADS"] = 12
println("Iniciando Script")
##Formato --> args = [experiments, iter, bal, prior, mov vecindario, tam vecindario, dmax, restrcprior, instance]
## Instance: PT = 0 ; PT_V1 = 1; PT_V2 = 2
## restrcprior: 0 significa que toma solo las prioridades de primer nivel. 1 significa que considera todos los niveles de prioridad.
## iter numero de iteraciones maximas
## el n/4 esta en duro. Cambiar!

########## PT_V1 ######################
#################################### PRIOR 10 BAL 0.5
#### 600 iteraciones ####

#args = [5, 600, 0.7, 20, [1 2 3], 3, 2500, 0, 0]
#include("./cultural/main.jl")

#args = [5, 600, 0.5, 20, [1 2 3], 3, 2500, 0, 0]
#include("./cultural/main.jl")

#args = [5, 600, 0.3, 20, [1 2 3], 3, 2500, 0, 0]
#include("./cultural/main.jl")
#dejar 
args = [5, 10, 1, 50, 1, 1, 0.1,  0.1, 0]
include("./cultural/main.jl")
#args = [5, 100, 10, 15, 1, 1, 0.05, 0.8, 0]
#include("./cultural/main.jl") 
# args = [1, 60, 6, 5, 1, 1, 0.05, 0.8, 0]
# include("./cultural/main.jl")
# args = [1, 80, 8, 5, 1, 1, 0.05, 0.8, 0]
# include("./cultural/main.jl")
# args = [1, 100, 10, 5, 1, 1, 0.05, 0.8, 0]
# include("./cultural/main.jl")
# args = [1, 120, 12, 5, 1, 1, 0.05, 0.8, 0]
# include("./cultural/main.jl")



# args = [1, 100, 10, 10, 1, 1, 0.05, 0.8, 0]
# include("./cultural/main.jl")


# args = [1, 10, 2, 5, 1, 1, 0.05, 0.8, 0]
# include("./cultural/main.jl")

# args = [1, 20, 2, 5, 1, 1, 0.05, 0.8, 0]
# include("./cultural/main.jl")

# args = [1, 20, 2, 5, 1, 1, 0.05, 0.8, 0]
# include("./cultural/main.jl")s

# args = [5, 20, 2, 5, 1, 1, 0.05, 0.8, 0]
# include("./cultural/main.jl")

#args = [5, 600, 0.7, 2, [1 2 3], 3, 2500, 0, 0] #Faltan estos
#include("./cultural/main.jl")

#args = [5, 600, 0.5, 2, [1 2 3], 3, 2500, 0, 0] # Faltan estos
#include("./cultural/main.jl")

#args = [5, 600, 0.3, 2, [1 2 3], 3, 2500, 0, 0] #Faltan estos
#include("./cultural/main.jl")

#################################################################

#################################################################
