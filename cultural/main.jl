global instancia = args[9]; #0 = PT, 1 = PT_V1, 2 = PT_V2
global last_obj = 0;

include("load_data.jl");
include("gurobi.jl");
include("helpers.jl");
include("cultural.jl");

using Statistics, DelimitedFiles;

global ide_exp = 3; # Numero del experimento.
global gi_order = []; # Orden de las grillas usadas

#Variables globales de balance y prioridad.
#Restricción de all prior # allprior = 0, allprior = 1;
#Criterio  # 1 = swap center random,  2 = swap center max distance,  3 = swap center prior bal

#Grilla
global M = clustering(); #grilla generada aleatoriamente
#global weight                 = load_grilla(ide_exp); #grilla cargada de la base de datos
# global M 		      = get_clusters(getindex.(findmax(weight,dims=2)[2],2))
# global M                = load_random_grilla(); #grilla cargada de manera aleatoria de la base de datos

M_aux = M;
#Matriz de adyacencia de zonas
#global adjacency_matrix = get_adjacency_matrix();

#grilla_db(); #generar una base de datos de grillas

#Matriz de conexiones
global c = connection_calculation();

#parámetros
balance = args[3];
prioridad = args[4];
dmax = args[7];
allprior = args[8];
tam_pob = 100
p_cross = 0.8  # Probabilidad de cruce
p_mut = 0.1  # Probabilidad de mutación
max_generaciones = 100  # Número máximo de generaciones
nr = 10  # Número de individuos a reemplazar

#= MAIN =#
# Variables metaheurística
    r_max = args[2]; #Número iteraciones.
neighborhood_structure = args[5]; #Tamaños estructuras de entorno.
len_N = args[6];  #Tamaño de los vecindarios.
k_max = length(neighborhood_structure); #k máximo.

#Numero de experimentos a realizar.
experimentos = args[1];

println("Parámetros");
println("   experimentos   ", args[1]);
println("   iteraciones    ", args[2]);
println("   balance        ", args[3]);
println("   prioridad      ", args[4]);
println("   mov_vecindario ", args[5]);
println("   tam_vecindario ", args[6]);
println("   dmax           ", args[7]);
println(" ");

#Limite de no mejoras.
const NO_IMPROVE_LIMIT = (r_max / 1);

objs_iter = 0;
objs_array = [];
exp_time_array = [];
println("algoritmo_cultural.")

println("Utilizando ", Threads.nthreads(), " hilo/s");


filename = "grilla_$ide_exp.txt";
mkpath("./grilla/")
open(joinpath("./grilla/", filename), "w") do file
    writedlm(file, M_aux)
end

for e = 1:experimentos
    C_test = zeros(Int64, length(CANDIDATAS))
    E_test = zeros(Int64, length(ESTACIONES))
    exp_time = @elapsed C_test, E_test, objs_iter, improve, obj = @time algoritmo_cultural(tam_pob, p_cross, p_mut, nr, NO_IMPROVE_LIMIT, e, max_generaciones)
    println("tiempo del experimento ", exp_time)
    append!(objs_array, objs_iter)
    append!(exp_time_array, exp_time)
    name = "$(balance)_$(prioridad)_exp_$(e)_$(r_max)_$(len_N)_$(improve)_$(obj)"
    filename = name * ".txt"
    open(joinpath("./grilla/$(nombre_instancia)/$(balance)_$(prioridad)/$(r_max)/$(neighborhood_structure)_$(len_N)/$(ide_exp)", filename), "a") do file
        write(file, "tiempo       = $exp_time \n")
    end

    filename = "matriz_conexiones_$(e).txt"
    open(joinpath("./Resultados finales/$(nombre_instancia)/$(balance)_$(prioridad)/$(r_max)/$(neighborhood_structure)_$(len_N)", filename), "a") do file
        #write(file, "c = $c_aux\n");
    end
    empty!(gi_order)
    global M = M_aux
end

let suma = 0.0, sumatimes = 0.0
    for i = 1:length(objs_array)
        suma = suma + objs_array[i]
    end
    for i = 1:length(exp_time_array)
        sumatimes = sumatimes + exp_time_array[i]
    end

    time_prom = sumatimes / length(exp_time_array)
    promedio = suma / length(objs_array)
    de = std(floor.(objs_array))
    best = minimum(objs_array)
    worst = maximum(objs_array)

    #Resumen resultados
    name = "result_exp_$(balance)_$(prioridad)_$(experimentos)_$(best)"
    filename = name * ".txt"
    open(joinpath("/home/guillermo/Desktop/Isaac/Resultados/fuzzy/$(nombre_instancia)/$(balance)_$(prioridad)/$(r_max)/$(neighborhood_structure)_$(len_N)/$(ide_exp)", filename), "w") do file
        write(file, "experimentos   = $experimentos \n")
        write(file, "promedio       = $promedio \n")
        write(file, "d.e            = $de   \n")
        write(file, "best           = $best \n")
        write(file, "worst          = $worst \n")
        write(file, "tiempo_prom    = $time_prom s\n")
        write(file, "tiempo_c/exp   = $exp_time_array\n")
        write(file, "n° grilla      = $ide_exp\n")
        write(file, "criterio       = swap center weight\n")
    end
end
