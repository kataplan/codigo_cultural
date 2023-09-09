
# Función para inicializar la población
function init_population(tam_pob)
    population = []
    for _ in 1:tam_pob
        individual = generar_individuo()  # Generar un individuo aleatorio
        push!(population, individual)
    end
    return population
end

function calculate_centroid(E, center)
    associated_stations = [i for i in 1:length(E) if E[i] == center]
    n = length(associated_stations)
    numerator = sum([dist[associated_stations[i], associated_stations[j]] for i in 1:n, j in 1:n])
    denominator = n * (n - 1) / 2
    centroid = numerator / denominator

    return centroid
end

function calculate_centroid_media(individual, estations, obj)
    centroid_sum = 0
    if (obj != Inf)
        for i in eachindex(individual)
            centroid = 0
            if individual[i] == 1
                centroid = calculate_centroid(estations, i)
            end
            centroid_sum += centroid
        end
    else

        return Inf
    end
    return centroid_sum / count(x -> x == 1, individual)

end

function fitness_population(not_fit_population)
    population = Vector{Any}(undef, length(not_fit_population))
    _obj_array =
        Threads.@threads for i in eachindex(not_fit_population)
            _obj, _E = evaluar_individuo(not_fit_population[i])
            population[i] = [not_fit_population[i], _E, _obj, calculate_centroid_media(not_fit_population[i], _E, _obj)]
        end
    println(population[1][4])

    return population
end


# Función para seleccionar padres mediante torneo binario
function seleccionar_padres(poblacion)
    padres = []
    while length(padres) < 2
        # Seleccionar dos individuos aleatorios de la población
        index_1 = rand(1:length(poblacion))
        index_2 = rand(1:length(poblacion))
        individuo1 = poblacion[index_1][1]
        individuo2 = poblacion[index_2][1]
        fitness1 = poblacion[index_1][3]
        fitness2 = poblacion[index_2][3]

        if fitness1 < fitness2
            push!(padres, individuo1)
        else
            push!(padres, individuo2)
        end
    end
    return padres
end

function selection(population, tam_pob)
    parents_pairs = []
    while length(parents_pairs) < tam_pob ÷ 2
        parents = seleccionar_padres(population)
        push!(parents_pairs, parents)
    end
    return parents_pairs
end


# Función para cruzar los padres y generar hijos
function crossover(pairs, crossover_tipe)
    children = []
    for pair in pairs
        padre1, padre2 = pair
        hijo1, hijo2 = realizar_cruce(padre1, padre2, crossover_tipe)
        push!(children, hijo1)
        push!(children, hijo2)
    end
    return children
end

function mutation(population, p_mut)
    mutated_population = copy(population)
    m = size(population, 1)
    n = size(population, 2)
    for inidividual in mutated_population
        for j in 1:n
            if rand() < p_mut
                # Realizar mutación en la característica j del individuo i
                mutated_population[i, j] = individual_mutation(inidividual, p_mut)
            end
        end
    end

    return mutated_population
end

function individual_mutation(individuo, prob_mutacion::Float64)
    nuevo_individuo = copy(individuo)
    # Verificar si el individuo es seleccionado para mutación
    for i in 1:length(individuo[1])
        # Verificar si el centro está designado en el individuo
        if individuo[1][i] == 1
            # Aplicar la probabilidad de mutación para cambiar el centro
            if rand() < prob_mutacion
                nuevo_individuo[1][i] = 0
                nuevo_individuo[1][rand(1:length(individuo))] = 1
            end
        end
    end
    return nuevo_individuo
end




# Función para generar un individuo aleatorio
function generar_individuo()
    C = init_solution_C_grid()
    return C
end

# Función para evaluar un individuo
function evaluar_individuo(individuo)
    return Gurobi_optimal(individuo)
end

# Función para realizar el cruce entre dos padres
function realizar_cruce(parent1, parent2, crossover_type)
    n = length(parent1)
    child1 = similar(parent1)
    child2 = similar(parent2)
    
    # Convertir los padres a máscaras
    maskParent1 = binary_to_mask(parent1)
    maskParent2 = binary_to_mask(parent2)

    if crossover_type == 1
        # Cruzamiento en un punto
        point = rand(2:length(maskParent1))
        child1[1:point] = maskParent1[1:point]
        child1[point+1:end] = maskParent2[point+1:end]
        child2[1:point] = maskParent2[1:point]
        child2[point+1:end] = maskParent1[point+1:end]
    elseif crossover_type == 2
        # Cruzamiento en dos puntos
        point1, point2 = sort(rand(2:length(maskParent1), 2))
        child1[1:point1] = maskParent1[1:point1]
        child1[point1+1:point2] = maskParent2[point1+1:point2]
        child1[point2+1:end] = maskParent1[point2+1:end]
        child2[1:point1] = maskParent2[1:point1]
        child2[point1+1:point2] = maskParent1[point1+1:point2]
        child2[point2+1:end] = maskParent2[point2+1:end]
    elseif crossover_type == 3
        # Cruzamiento uniforme
        for i in 1:n
            if rand() < 0.5
                child1[i] = maskParent1[i]
                child2[i] = maskParent2[i]
            else
                child1[i] = maskParent2[i]
                child2[i] = maskParent1[i]
            end
        end
    else
        error("Tipo de cruzamiento no válido")
    end

    return child1, child2
end

function binary_to_mask(binary)
    print(binary)
    # Obtener las posiciones de los "1" en el vector binario
    mask = findall(x -> x == 1, binary)
    print(mask)
    return mask
end


function acceptance(belief_network, population, max_size_belefief_space)
    lenght_belief_space = size(belief_network, 1)
    for individual in population
        for j in 1:lenght_belief_space
            if ((individual[j+2] >= belief_network[j, 2] && individual[j+2] <= belief_network[j, 3]) || individual[j+2] <= belief_network[j, 2])
                println("Inidividuo agregado al espacio de creencias. ")
                # El individuo cumple los requisitos del componente normativo
                if length(belief_network[j, 1]) >= max_size_belefief_space
                    # Si la memoria cultural está llena, eliminamos el último individuo
                    belief_network[j, 1] = belief_network[j, 1][1:end-1]
                end

                # Agregamos el nuevo individuo a la memoria cultural
                push!(belief_network[j, 1], individual)
                # Ordenamos el array belief_network[j, 1] en base al fitness de los individuos (mayor a menor)
                sort!(belief_network[j, 1], by=x -> x[j+2], rev=false)

                # Actualizamos los rangos del intervalo I y los puntajes L y U
                belief_network[j, 2] = belief_network[j, 1][1][j+2]
                belief_network[j, 3] = belief_network[j, 1][end][j+2]
            end
        end
    end

    return belief_network
end


function init_belief_network(n)

    belief_network = Array{Any}(undef, n, 3)

    for i in 1:n
        belief_network[i, 1] = []
        belief_network[i, 2] = -Inf
        belief_network[i, 3] = Inf
    end

    return belief_network
end




function get_best(population)
    best_fitness = minimum([individual[3] for individual in population])  # Obtener el fitness máximo en la columna 2
    best_individual_index = findfirst(x -> x[3] == best_fitness, population)  # Encontrar el índice del fitness máximo
    best = population[best_individual_index]
    return best  # Devolver el individuo con el fitness máximo
end

function influence(population, belief_network)
    new_population = population  # Inicializar la nueva población con la población original
    n = size(belief_network, 1)
    for i in 1:n
        individuals = belief_network[i, 1]  # Obtener los individuos de belief_network[i, 1]
        new_population = vcat(new_population, individuals)  # Concatenar horizontalmente a la nueva población
    end
    return new_population
end

# Algoritmo genético con memoria cultural
function algoritmo_cultural(tam_pob, p_cross, p_mut, max_generaciones, max_size_belefief_space, crossover_tipe)

    no_fit_population = init_population(tam_pob)
    population = fitness_population(no_fit_population)
    belief_network = init_belief_network(length(population[1][3:end]))
    belief_network = acceptance(belief_network, population, max_size_belefief_space)
    best_individual = get_best(population)
    i = 0
    best_generation = 0
    while i < max_generaciones
        println("-----Actual Generacion ", i, " -----")
        population = influence(population, belief_network)
        p = selection(population, tam_pob)
        ti = crossover(p, crossover_tipe)
        ti = mutation(ti, p_mut)
        population = fitness_population(ti)
        belief_network = acceptance(belief_network, population, max_size_belefief_space)
        best_individual_i = get_best(population)
        if best_individual[3] > best_individual_i[3]
            best_individual = best_individual_i
            best_generation = i
        end
        i += 1
    end
    return best_individual, best_generation
end

