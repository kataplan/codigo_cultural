global instancia = args[9]; #0 = PT, 1 = PT_V1, 2 = PT_V2
global last_obj = 0;

include("load_data.jl");
include("gurobi.jl");
include("helpers.jl");
include("cultural.jl");
include("cultural_helpers.jl");

using Statistics, DelimitedFiles;

global ide_exp = 3; # Numero del experimento.
global gi_order = []; # Orden de las grillas usadas

#Variables globales de balance y prioridad.
#Restricción de all prior # allprior = 0, allprior = 1;
#Criterio  # 1 = swap center random,  2 = swap center max distance,  3 = swap center prior bal

#Grilla
global M = clustering(); #grilla generada con kmenas Size 15x452
# global weight                 = load_grilla(ide_exp); #grilla cargada de la base de datos
# global M 		      = get_clusters(getindex.(findmax(weight,dims=2)[2],2))

# global M                = load_random_grilla(); #grilla cargada de manera aleatoria de la base de datos

M_aux = M;
#Matriz de adyacencia de zonas
#global adjacency_matrix = get_adjacency_matrix();

#grilla_db(); #generar una base de datos de grillas

#Matriz de conexiones
global c = connection_calculation();

#Paramaters
experimentos = args[1]
cross_size = args[2]
max_size_belefief_space = args[3]  # Número de individuos a reemplazar
max_generations = args[4]  # Número máximo de generaciones
crossover_type = args[5]
mutation_size = args[6]  # Probabilidad de mutación
influence_size = args[7]  # Probabilidad de cruce
end_rule = args[8]

println("Parámetros");
println("   Tamaño de poblacion           = ", args[2]);
println("   Indiviudos destacados máximos = ", args[3]);
println("   Generaciones Máximas          = ", args[4]);
println("   Tipo del crossover            = ", args[5]);
println("   Individuos para mutación      = ", args[6]);
println("   Individuos para influencia    = ", args[7]);
println("   Regla de termino              = ", args[8]);
println(" ");

#Limite de no mejoras.
#const NO_IMPROVE_LIMIT = (r_max / 1);

objs_array = [];
exp_time_array = [];
println("algoritmo_cultural.")
ENV["JULIA_NUM_THREADS"] = 12
println("Hilos disponibles ", Threads.nthreads(), " hilo/s");

# filename = "grilla_$ide_exp.txt";
# mkpath("./grilla/")
# open(joinpath("./grilla/", filename), "w") do file
#     writedlm(file, M_aux)
# end

for e = 1:experimentos
    C_test = zeros(Int64, length(CANDIDATAS))
    E_test = zeros(Int64, length(ESTACIONES))
    exp_time = @elapsed individuo, generacion, best_individuals_per_generation, best_overall_individuals = @time cultural_algorithm(cross_size, influence_size, mutation_size, max_generations, max_size_belefief_space, crossover_type, end_rule)
    println("tiempo del experimento : ", exp_time)
    println("individuo              : ", individuo["individual"])
    println("E                      : ", individuo["E"])
    println("fitness                : ", individuo["obj"])
    println("best generation        : ", generacion)

    c_aux = individuo["individual"]
    append!(objs_array, individuo["obj"])
    append!(exp_time_array, exp_time)
    name = "$(cross_size)_$(max_generations)_$(crossover_type)_$(mutation_size)_$(max_size_belefief_space)_$(e)"
    filename = name * ".txt"
    open(joinpath("cultural/resultados", filename), "w") do file
        write(file,"tiempo       = $(exp_time) \nfitness                : , $(individuo["obj"])\nBest overall:\n$(best_overall_individuals)\nBest per generation:\n$(best_individuals_per_generation)")
    end

    # filename = "matriz_conexiones_$(e).txt"
    # open(joinpath("./Resultados finales/$(tam_pob)_$(max_generations)_$(crossover_tipe)_$(mutation_tipe)_$(max_size_belefief_space)", filename), "w") do file
    #     write(file, "c = $c_aux\n");
    # end
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
    open(joinpath("cultural/resultados_text", filename), "w") do file
        write(file, "experimentos      = $experimentos \n")
        write(file, "promedio          = $promedio \n")
        write(file, "d.e               = $de   \n")
        write(file, "best              = $best \n")
        write(file, "worst             = $worst \n")
        write(file, "tiempo_prom       = $time_prom s\n")
        write(file, "tiempo_c/exp      = $exp_time_array\n")
        write(file, "n° grilla         = $ide_exp\n")
        write(file, "end number        = $max_generations\n")
        write(file, "end condition     = $end_rule\n")
        write(file, "beleif space size = $max_size_belefief_space\n")
        write(file, "population size   = $cross_size\n")
        write(file, "influence size    = $influence_size\n")
        write(file, "mutation size     = $mutation_size\n")
        write(file, "crossover type    = $crossover_type\n")
    end
end
