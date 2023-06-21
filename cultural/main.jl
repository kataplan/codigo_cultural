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
#global M 		      = get_clusters(getindex.(findmax(weight,dims=2)[2],2))

# global M                = load_random_grilla(); #grilla cargada de manera aleatoria de la base de datos

M_aux = M;
#Matriz de adyacencia de zonas
#global adjacency_matrix = get_adjacency_matrix();

#grilla_db(); #generar una base de datos de grillas

#Matriz de conexiones
global c = connection_calculation();

#Paramaters
experimentos = args[1]
tam_pob = args[2]
max_size_belefief_space = args[3]  # Número de individuos a reemplazar
max_generaciones = args[4]  # Número máximo de generaciones
crossover_tipe = args[5]
mutation_tipe = args[6]
p_mut = args[7]  # Probabilidad de mutación
p_cross = args[8]  # Probabilidad de cruce

println("Parámetros");
println("   Tamaño de poblacion           = ", args[2]);
println("   Indiviudos destacados máximos = ", args[3]);
println("   Generaciones Máximas          = ", args[4]);
println("   Tipo del crossover            = ", args[5]);
println("   Tipo de la mutación           = ", args[6]);
println("   Probabilidad de mutación      = ", args[7]);
println("   Probabilidad de cruce         = ", args[8]);
println(" ");

#Limite de no mejoras.
#const NO_IMPROVE_LIMIT = (r_max / 1);

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
    exp_time = @elapsed individuo, generacion = @time algoritmo_cultural(tam_pob, p_cross, p_mut, max_generaciones, max_size_belefief_space, crossover_tipe)
    println("tiempo del experimento : ", exp_time)
    #println("individuo              : ", individuo[1])
    println("fitness                : ", individuo[2])
    #println("E                      : ", individuo[3])
    println("generación             : ", generacion)
    append!(objs_array, objs_iter)
    append!(exp_time_array, exp_time)
    name = "$(tam_pob)_$(max_generaciones)_crossover_$(crossover_tipe)_$(mutation_tipe)_$(p_cross)_$(improve)_$(obj)"
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
