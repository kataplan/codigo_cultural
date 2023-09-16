# Function to initialize the population
function init_population(pop_size)
    population = []
    for _ in 1:pop_size
        individual = generate_individual()  # Generate a random individual
        push!(population, individual)
    end
    return population
end

function calculate_centroid(E, center)
    associated_stations = findall(x -> x == center, E)
    n = length(associated_stations)
    numerator = sum([dist[associated_stations[i], associated_stations[j]] for i in 1:n, j in 1:n])
    denominator = n * (n - 1) / 2
    centroid = numerator / denominator

    return centroid
end

function calculate_centroid_mean(individual, stations, obj)
    centroid_sum = 0
    if (obj != Inf)
        for i in eachindex(individual)
            centroid = 0
            if individual[i] == 1
                centroid = calculate_centroid(stations, i)
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
    Threads.@threads for i in eachindex(not_fit_population)
        _obj, _E = evaluate_individual(not_fit_population[i])
        population[i] = [not_fit_population[i], _E, _obj, calculate_centroid_mean(not_fit_population[i], _E, _obj)]
    end
    return population
end

# Function to select parents using binary tournament
function select_parents(population)
    parents = []
    while length(parents) < 2
        # Select two random individuals from the population
        index_1 = rand(1:length(population))
        index_2 = rand(1:length(population))
        individual1 = population[index_1][1]
        individual2 = population[index_2][1]
        fitness1 = population[index_1][3]
        fitness2 = population[index_2][3]

        if fitness1 < fitness2
            push!(parents, individual1)
        else
            push!(parents, individual2)
        end
    end
    return parents
end

function selection(population, pop_size)
    parents_pairs = []
    while length(parents_pairs) < pop_size ÷ 2
        parents = select_parents(population)
        push!(parents_pairs, parents)
    end
    return parents_pairs
end

# Function to cross parents and generate children
function crossover(pairs, crossover_type)
    println("crossover Process starting")
    children = []
    for pair in pairs
        parent1, parent2 = pair
        child1, child2 = perform_crossover(parent1, parent2, crossover_type)
        push!(children, child2)
        push!(children, child1)

    end
    return children
end

function mutation(population, p_mut)
    println("Mutation Process starting")
    mutated_population = copy(population)
    m = size(population, 1)
    n = size(population, 2)
    for individual in mutated_population
        if (rand(0:1) < p_mut)
            individual_mutation(individual)
        end
    end
    return mutated_population
end

function select_random_position(arr, value::Int)
    valid_positions = findall(x -> x == 1 && x != value, arr)
    if isempty(valid_positions)
        return nothing  # No se encontraron posiciones válidas
    end

    random_position = valid_positions[rand(1:length(valid_positions))]
    return random_position
end

function individual_mutation(individual)
    new_individual = copy(individual)
    ones = count(x -> x == 1, new_individual)
    if (ones != 15)
        error("error ENTRO MAL numero de unos>", ones, "      ", new_individual)
    end
    # Check if the individual is selected for mutation
    random_cluster = rand(1:15)
    # Get the indices of centers in the same cluster
    cluster_centers = M[random_cluster, :]
    matching_positions = 100000
    for i in 1:length(cluster_centers)
        if cluster_centers[i] == 1 && new_individual[i] == 1
            matching_positions = i
            break
        end
    end
    new_individual[matching_positions] = 0
    mutation_position = select_random_position(cluster_centers, matching_positions)

    new_individual[mutation_position] = 1
    ones = count(x -> x == 1, new_individual)

    if (ones != 15)
        error("error MALA MUTACION numero de unos>", ones, "      ", new_individual)
    end
    return new_individual
end



function check_ones_count(individual)
    if count(x -> x == 1, individual) != 15
        return false
    end
    return true
end

# Function to generate a random individual
function generate_individual()
    C = init_solution_C_grid()
    return C
end

# Function to evaluate an individual
function evaluate_individual(individual)
    return Gurobi_optimal(individual)
end

# Function to perform crossover between two parents
function perform_crossover(parent1, parent2, crossover_type)
    n = length(parent1)
    # Convert parents to masks
    maskParent1 = binary_to_mask(parent1, M)
    maskParent2 = binary_to_mask(parent2, M)

    child1 = similar(maskParent1)
    child2 = similar(maskParent2)
    if crossover_type == 1
        # Single-point crossover
        point = rand(2:length(maskParent1))
        child1[1:point] = maskParent1[1:point]
        child1[point+1:end] = maskParent2[point+1:end]
        child2[1:point] = maskParent2[1:point]
        child2[point+1:end] = maskParent1[point+1:end]
    elseif crossover_type == 2
        # Two-point crossover
        point1, point2 = sort(rand(2:length(maskParent1), 2))
        child1[1:point1] = maskParent1[1:point1]
        child1[point1+1:point2] = maskParent2[point1+1:point2]
        child1[point2+1:end] = maskParent1[point2+1:end]
        child2[1:point1] = maskParent2[1:point1]
        child2[point1+1:point2] = maskParent1[point1+1:point2]
        child2[point2+1:end] = maskParent2[point2+1:end]
    elseif crossover_type == 3
        # Uniform crossover
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
        error("Invalid crossover type")
    end

    childUnmask1 = mask_to_binary(child1, n)
    childUnmask2 = mask_to_binary(child2, n)

    return childUnmask1, childUnmask2
end

function mask_to_binary(child_mask, parent_size)
    binary_child = zeros(Int, parent_size)
    for pos in child_mask
        binary_child[pos] = 1
    end
    return binary_child
end

function binary_to_mask(binary, M)
    # Obtener las posiciones de los "1" en el vector binario
    mask = findall(x -> x == 1, binary)
    # Ordenar las posiciones en función de la fila correspondiente en M
    sorted_mask = sort(mask, by=x -> findfirst(isequal(1), M[:, x]))
    return sorted_mask
end

function acceptance(belief_network, population, max_belief_space_size)
    belief_space_length = size(belief_network, 1)
    for individual in population
        for j in 1:belief_space_length
            if ((individual[j+2] > belief_network[j, 2] && individual[j+2] < belief_network[j, 3]) || individual[j+2] < belief_network[j, 2])
                
                println("Individual added to the belief space.")
                println(j)
                println(belief_network[j, 2])
                println(individual[j+2])
                println(belief_network[j, 3])
                # The individual meets the normative component requirements
                if length(belief_network[j, 1]) >= max_belief_space_size
                    # If the cultural memory is full, remove the last individual
                    belief_network[j, 1] = belief_network[j, 1][1:end-1]
                end

                # Add the new individual to the cultural memory
                push!(belief_network[j, 1], individual)
                # Sort the belief_network[j, 1] array based on individual fitness (highest to lowest)
                sort!(belief_network[j, 1], by=x -> x[j+2], rev=false)

                # Update the range of interval I and the L and U scores
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
    best_fitness = minimum([individual[3] for individual in population])  # Get the maximum fitness in column 2
    best_individual_index = findfirst(x -> x[3] == best_fitness, population)  # Find the index of the maximum fitness
    best = population[best_individual_index]
    return best  # Return the individual with the maximum fitness
end

function influence(population, belief_network)
    new_population = population  # Initialize the new population with the original population
    n = size(belief_network, 1)
    for i in 1:n
        individuals = belief_network[i, 1]  # Get the individuals from belief_network[i, 1]
        new_population = vcat(new_population, individuals)  # Concatenate horizontally to the new population
    end
    return new_population
end



# Cultural algorithm
function cultural_algorithm(pop_size, p_cross, p_mut, max_generations, max_belief_space_size, crossover_type)
    ENV["JULIA_NUM_THREADS"] = 12
    non_fit_population = init_population(pop_size)
    population = fitness_population(non_fit_population)
    belief_network = init_belief_network(length(population[1][3:end]))
    belief_network = acceptance(belief_network, population, max_belief_space_size)
    best_individual = get_best(population)
    i = 0
    best_generation = 0
    while i < max_generations - 1
        println("-----Current Generation ", i + 1, " -----")
        if(max_generations*0.1 < i)
            population = influence(population, belief_network)
        end
        p = selection(population, pop_size)
        ti = crossover(p, crossover_type)
        check_ones_count(ti)
        ti = mutation(ti, p_mut)
        check_ones_count(ti)
        population = fitness_population(ti)
        belief_network = acceptance(belief_network, population, max_belief_space_size)
        best_individual_i = get_best(population)
        if best_individual[3] > best_individual_i[3]
            println("new Best individual in ", i+1, " with: ", best_individual[3])
            best_individual = best_individual_i
            best_generation = i
        end
        i += 1
    end
    return best_individual, best_generation
end
