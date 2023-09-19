# Function to initialize the population
condition_array = ["obj", "centroid"]

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
        _obj, _E, slack_array = evaluate_individual(not_fit_population[i])
        population[i] = Dict("individual" => not_fit_population[i],
            "E" => _E,
            "obj" => _obj,
            "centroid" => calculate_centroid_mean(not_fit_population[i], _E, _obj),
            "slack" => slack_array)
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
        if population[index_1]["obj"] < population[index_2]["obj"]
            push!(parents, population[index_1]["individual"])
        else
            push!(parents, population[index_2]["individual"])
        end
    end

    return parents
end

function selection(population, pop_size, influence_percentage, mutation_percentage)
    # Calcular el tamaño de la parte restante
    influence_size = Int(ceil(pop_size * influence_percentage))
    mutation_size = Int(ceil(pop_size * mutation_percentage))
    remaining_size = Int(pop_size - influence_size - mutation_size)

    # Crear índices desordenados para la población
    indices = randperm(Int(pop_size))
    println(length(population))
    # Dividir la población en las tres partes
    influence_population = [population[i] for i in indices[1:influence_size]]
    mutation_population = [population[i] for i in indices[influence_size+1:influence_size+mutation_size]]
    remaining_population = [population[i] for i in indices[influence_size+mutation_size+1:Int(pop_size)]]
   

    return remaining_population, influence_population, mutation_population
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

function mutation(population)
    println("Mutation Process starting")
    remaining_pairs = []

  
    mutated_population = []
    for individual in population
        mutated_population = push!(mutated_population, individual_mutation(individual))
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
    new_individual = copy(individual["individual"])

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
    if (mutation_position !== nothing)
        new_individual[mutation_position] = 1
    end


    return new_individual
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
    maskParent1 = binary_to_mask(parent1)
    maskParent2 = binary_to_mask(parent2)
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

function binary_to_mask(binary)
    # Obtener las posiciones de los "1" en el vector binario
    mask = findall(x -> x == 1, binary)
    # Ordenar las posiciones en función de la fila correspondiente en M
    sorted_mask = sort(mask, by=x -> findfirst(isequal(1), M[:, x]))
    return sorted_mask
end

function acceptance(belief_network, population, max_belief_space_size)
    for individual in population
        for condition in condition_array
            if ((individual[condition] > belief_network[condition]["lower_bound"] && individual[condition] < belief_network[condition]["upper_bound"]) || individual[condition] < belief_network[condition]["lower_bound"])
                println("Individual added to the belief space.")
                # The individual meets the normative component requirements
                if length(belief_network[condition]["individuals"]) >= max_belief_space_size
                    # If the cultural memory is full, remove the last individual
                    belief_network[condition]["individuals"] = belief_network[condition]["individuals"][1:end-1]
                end

                # Add the new individual to the cultural memory
                push!(belief_network[condition]["individuals"], individual)
                # Sort the belief_network[j, 1] array based on individual fitness (highest to lowest)
                sort!(belief_network[condition]["individuals"], by=x -> x[condition], rev=false)

                # Update the range of interval I and the L and U scores
                belief_network[condition]["lower_bound"] = belief_network[condition]["individuals"][1][condition]
                belief_network[condition]["upper_bound"] = belief_network[condition]["individuals"][end][condition]
            end
        end
    end

    return belief_network
end

function init_belief_network()

    belief_network = Dict{Any,Dict{String,Any}}()  # Define the belief_network as an empty dictionary

    for condition in condition_array
        belief_network[condition] = Dict("individuals" => [], "lower_bound" => -Inf, "upper_bound" => Inf)
    end

    return belief_network
end

function get_best(population)
    best_fitness = minimum([individual["obj"] for individual in population])  # Get the maximum fitness in column 2
    best_individual_index = findfirst(x -> x["obj"] == best_fitness, population)  # Find the index of the maximum fitness
    best = population[best_individual_index]
    return best  # Return the individual with the maximum fitness
end

function influence(population, belief_network, influence_percentage)
    modified_individuals = []

    # Calculate the number of individuals to select (10% of the population)
    num_to_select = Int(ceil(influence_percentage * length(population)))
    # Randomly select num_to_select individuals from the current population
    selected_individuals = sample(population, num_to_select, replace=false)
    condition_index = 1
    # Iterate over the selected individuals
    for individual in selected_individuals
        # Iterate over the conditions in the belief network
        influencer = belief_network[condition_array[condition_index]]["individuals"][rand(1:length(belief_network[condition_array[condition_index]]["individuals"]))]
        mask_influencer = binary_to_mask(influencer["individual"])
        mask_individual = binary_to_mask(individual["individual"])
        for i in 1:length(individual["slack"])
            if (individual["slack"][i] < influencer["slack"][i])
                mask_individual[i] = mask_influencer[i]
            end
        end
        modified_individuals = push!(modified_individuals, mask_to_binary(mask_individual, length(individual["individual"])))
        # Increment the slack_index, cycling through condition_array
        condition_index = mod(condition_index, length(condition_array)) + 1
    end

    return modified_individuals
end

# Helper function to calculate slack for each cluster
function calculate_slack(slack_array, M)


    return slack_values
end


# Cultural algorithm
function cultural_algorithm(pop_size, p_influ, p_mut, max_generations, max_belief_space_size, crossover_type, experiment)
    ENV["JULIA_NUM_THREADS"] = 12
    non_fit_population = init_population(pop_size)
    population = fitness_population(non_fit_population)
    belief_network = init_belief_network()
    belief_network = acceptance(belief_network, population, max_belief_space_size)
    best_individual = get_best(population)
    i = 0
    best_generation = 0
    while i < max_generations - 1
        println("-----Current Generation ", i + 1, " -----")
        #if (max_generations * 0.1 < i)
        #    population = influence(population, belief_network, p_cross)
        #end
        ti_cross, ti_mut, ti_influence = selection(population, pop_size, p_influ, p_mut)
        ti_cross = crossover(ti_cross, crossover_type)
        ti_influence = influence(ti_influence, belief_network, p_influ)
        ti_mut = mutation(ti_mut)
        println(length(ti_cross))
        println(length(ti_influence))
        println(length(ti_mut))
        population = fitness_population(vcat(ti_cross, ti_mut, ti_influence))
        belief_network = acceptance(belief_network, population, max_belief_space_size)
        best_individual_i = get_best(population)
        if best_individual["obj"] > best_individual_i["obj"]
            println("new Best individual in ", i + 1, " with: ", best_individual["obj"])
            best_individual = best_individual_i
            best_generation = i
        end
        i += 1
    end
    name = "result_exp_$(experiment)_resultados"
    filename = name * ".txt"
    open(joinpath("cultural/resultados_text", filename), "w") do file
        write(file, "new Best individual in   = $best_generation \n")
        write(file, "Best individual          = $best_individual \n")

    end
    return best_individual, best_generation
end
