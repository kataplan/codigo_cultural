# Función para inicializar la población
function init_population(tam_pob)
    population = []
    for _ in 1:tam_pob
        individual = generar_individuo()  # Generar un individuo aleatorio
        push!(population, individual)
    end
    return population
end

function fintess_population(not_fit_population)
    population = []
    for no_fit_individual in not_fit_population
        _obj, _E = evaluar_individuo(no_fit_individual)
        individual = [no_fit_individual, _obj, _E]
        push!(population, individual)
    end
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
        fitness1 = poblacion[index_1][2]
        fitness2 = poblacion[index_2][2]
        
        if fitness1 > fitness2
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
    for i in 1:m
        for j in 1:n
            if rand() < p_mut
                # Realizar mutación en la característica j del individuo i
                mutated_population[i, j] = realizar_mutacion(population[i, j])
            end
        end
    end

    return mutated_population
end

# Función para reemplazar los individuos menos aptos
function reemplazar_individuos(poblacion, hijos, nr)
    ordenar_poblacion(poblacion)
    ordenar_poblacion(hijos)
    for i in 1:nr
        poblacion[end-i+1] = hijos[i]
    end
end



# Función para evaluar el criterio de parada
function criterio_de_parada_cumplido(mejor_solucion, NO_IMPROVE_LIMIT)
    # Implementar el criterio de parada adecuado a tu problema
    if mejor_solucion[:aptitud] >= NO_IMPROVE_LIMIT
        return true
    else
        return false
    end
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

    if crossover_type == 1
        # Cruzamiento en un punto
        point = rand(2:n-1)
        child1[1:point] = parent1[1:point]
        child1[point+1:end] = parent2[point+1:end]
        child2[1:point] = parent2[1:point]
        child2[point+1:end] = parent1[point+1:end]
    elseif crossover_type == 2
        # Cruzamiento en dos puntos
        point1, point2 = sort(rand(2:n-1, 2))
        child1[1:point1] = parent1[1:point1]
        child1[point1+1:point2] = parent2[point1+1:point2]
        child1[point2+1:end] = parent1[point2+1:end]
        child2[1:point1] = parent2[1:point1]
        child2[point1+1:point2] = parent1[point1+1:point2]
        child2[point2+1:end] = parent2[point2+1:end]
    elseif crossover_type == 3
        # Cruzamiento uniforme
        for i in 1:n
            if rand() < 0.5
                child1[i] = parent1[i]
                child2[i] = parent2[i]
            else
                child1[i] = parent2[i]
                child2[i] = parent1[i]
            end
        end
    else
        error("Tipo de cruzamiento no válido")
    end

    return child1, child2
end



function acceptance(belief_network, population, max_size_belefief_space)
    lenght_belief_space = size(belief_network, 1)
    for individual in population
        for j in 1:lenght_belief_space
            if ((individual[2] >= belief_network[j, 2] && individual[2] <= belief_network[j, 3]) || individual[2] >= belief_network[j, 3])
                println("Inidividuo agregado al espacio de creencias. ")
                # El individuo cumple los requisitos del componente normativo
                if length(belief_network[j, 1]) >= max_size_belefief_space
                    # Si la memoria cultural está llena, eliminamos el último individuo
                    belief_network[j, 1] = belief_network[j, 1][1:end-1]
                end

                # Agregamos el nuevo individuo a la memoria cultural
                push!(belief_network[j, 1], individual)
                # Ordenamos el array belief_network[j, 1] en base al fitness de los individuos (mayor a menor)
                sort!(belief_network[j, 1], by=x -> x[2], rev=true)

                # Actualizamos los rangos del intervalo I y los puntajes L y U
                belief_network[j, 2] = individual[2]
                belief_network[j, 3] = individual[2]
            end
        end
    end
    return belief_network
end


function init_belief_network(max_size_belefief_space)
    #la dimension es uno xq solo estoy usando el fitness
    n = 1
    belief_network = Array{Any}(undef, n, 3)

    for i in 1:n
        belief_network[i, 1] = []
        belief_network[i, 2] = -Inf
        belief_network[i, 3] = Inf
    end

    return belief_network
end

function isinside(intervalo, valor)
    min_val, max_val = intervalo
    return min_val <= valor <= max_val
end

function update_dimension(intervalo, L, U, individuo)
    min_val, max_val = intervalo
    if individuo < min_val
        min_val = individuo
        L = performance(individuo)  # Calcula el puntaje de performance del valor mínimo
    end
    if individuo > max_val
        max_val = individuo
        U = performance(individuo)  # Calcula el puntaje de performance del valor máximo
    end
end
# Función para explorar culturalmente
function explorar_culturalmente(memoria_cultural, e)
    exploradores = memoria_cultural[rand(1:length(memoria_cultural), e)]
    for explorador in exploradores
        realizar_exploracion(explorador)
    end
end

# Función para imprimir información de la generación
function imprimir_informacion_generacion(gen, mejor_solucion)
    println("Generación: ", gen)
    println("Mejor solución encontrada: ", mejor_solucion)
    println("--------------------")
end

function get_maximum(population)
    best_fitness = maximum([individual[2] for individual in population])  # Obtener el fitness máximo en la columna 2
    best_individual_index = findfirst(x -> x[2] == best_fitness, population)  # Encontrar el índice del fitness máximo
    return population[best_individual_index]  # Devolver el individuo con el fitness máximo
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
    population = fintess_population(no_fit_population)
    belief_network = init_belief_network(max_size_belefief_space)
    belief_network = acceptance(belief_network, population, max_size_belefief_space)
    best_individual = get_maximum(population)
    i = 0

    while i < max_generaciones
        population = influence(population, belief_network)
        p = selection(population, tam_pob)
        ti = crossover(p, crossover_tipe)
        #ti = mutation(ti, p_mut)
        population = fintess_population(ti)
        belief_network = acceptance(belief_network, population, max_size_belefief_space)
        best_individual_i = get_maximum(population)
        if best_individual[2] > best_individual_i[2]
            best_individual = best_individual_i
        end
        i += 1
    end

    return best_individual[3]
end

