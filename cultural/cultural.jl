# Function to initialize the population
condition_array = ["obj"]

function init_population(pop_size)
    population = []
    for _ in 1:pop_size
        individual = init_solution_C_grid()  # Generate a random individual
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
    centroid_sum = 0.0
    station_indices = findall(individual .== 1)
    num_selected_stations = length(station_indices)

    if obj == Inf
        return Inf
    end

    if num_selected_stations > 0
        # Divide los índices de las estaciones en minibatches de un tamaño fijo, por ejemplo, 100 estaciones por minibatch.
        minibatch_size = 50
        num_minibatches = ceil(Int, num_selected_stations / minibatch_size)

        for i = 1:num_minibatches
            start_idx = (i - 1) * minibatch_size + 1
            end_idx = min(i * minibatch_size, num_selected_stations)
            minibatch_indices = station_indices[start_idx:end_idx]

            if !isempty(minibatch_indices)
                minibatch_centroids = [calculate_centroid(stations, center) for center in minibatch_indices]
                centroid_sum += sum(minibatch_centroids)
            end
        end
    end

    return centroid_sum / num_selected_stations
end


function fitness_population(not_fit_population, count_matrix)
    population = Vector{Any}(undef, length(not_fit_population))
    Threads.@threads for i in eachindex(not_fit_population)
        _obj, _E, slack_array = Gurobi_optimal(not_fit_population[i])
        population[i] = Dict(
            "individual" => not_fit_population[i],
            "E" => _E,
            "obj" => _obj,
            "centroid" => calculate_centroid_mean(not_fit_population[i], _E, _obj),
            "slack" => slack_array
        )
        count_matrix = update_center_counts(not_fit_population[i], count_matrix)
    end
    return population, count_matrix
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

function selection(population, population_sizes)
    # Crear índices desordenados para la población
    indices = randperm(population_sizes["pop_size"])
    # Dividir la población en las tres partes
    influence_population = [population[i] for i in indices[1:population_sizes["influence"]]]
    mutation_population = [population[i] for i in indices[population_sizes["influence"]+1:population_sizes["influence"]+population_sizes["mutation"]]]
    remaining_population = [population[i] for i in indices[population_sizes["influence"]+population_sizes["mutation"]+1:population_sizes["pop_size"]]]
    return remaining_population, influence_population, mutation_population
end


# Function to cross parents and generate children
function crossover(population, crossover_type)
    println("Crossover Process Starting")
    children = []
    pairs = []
    while length(pairs) < length(population) ÷ 2
        parents = select_parents(population)
        push!(pairs, parents)
    end
    for pair in pairs
        parent1, parent2 = pair
        child1, child2 = perform_crossover(parent1, parent2, crossover_type)
        push!(children, child2)
        push!(children, child1)
    end
    return children
end

function mutation(population, center_counts)
    println("Mutation Process starting")
    mutated_population = []
    for individual in population
        mutated_population = push!(mutated_population, individual_mutation(individual, center_counts))
    end
    return mutated_population
end


function individual_mutation(individual, center_counts)
    new_individual = copy(individual["individual"])
    mask_individual = binary_to_mask(new_individual)
    count_array = calculate_count_array(mask_individual, center_counts)
    # Encontrar la variedad máxima y los índices de los centros con la misma variedad
    max_count = maximum(count_array)
    max_count_indices = findall(count_array .== max_count)
    # Elegir un centro al azar entre los que tienen la máxima variedad
    if !isempty(max_count_indices)
        random_center_index = rand(max_count_indices)
        mutation_position = select_random_position(M[random_center_index, :], mask_individual[random_center_index])
        if mutation_position !== nothing
            mask_individual[random_center_index] = mutation_position
        end
    end
    new_individual = mask_to_binary(mask_individual, length(individual["individual"]))
    return new_individual
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
        belief_network[condition] =
            Dict(
                "individuals" => [],
                "lower_bound" => -Inf,
                "upper_bound" => Inf
            )
    end

    return belief_network
end

function get_best(population)
    best_fitness = minimum([individual["obj"] for individual in population])  # Get the maximum fitness in column 2
    best_individual_index = findfirst(x -> x["obj"] == best_fitness, population)  # Find the index of the maximum fitness
    best = population[best_individual_index]
    return best  # Return the individual with the maximum fitness
end

function influence(population, belief_network)
    println("Influence Process starting")
    modified_individuals = []
    condition_index = 1
    # Iterate over the selected individuals
    for individual in population
        # Iterate over the conditions in the belief network
        influencer = belief_network[condition_array[condition_index]]["individuals"][rand(1:length(belief_network[condition_array[condition_index]]["individuals"]))]
        mask_influencer = binary_to_mask(influencer["individual"])
        mask_individual = binary_to_mask(individual["individual"])
        for i in 1:length(individual["slack"])
            if (individual["slack"][i] < influencer["slack"][i])
                mask_individual[i] = mask_influencer[i]
            end
        end
        new_individual = mask_to_binary(mask_individual, length(individual["individual"]))
        #end debug
        modified_individuals = push!(modified_individuals, new_individual)
        # Increment the slack_index, cycling through condition_array
        condition_index = mod(condition_index, length(condition_array)) + 1
    end

    return modified_individuals
end


# Cultural algorithm
function cultural_algorithm(cross_size, influence_size, mutation_size, end_number, max_belief_space_size, crossover_type, end_rule)
    ENV["JULIA_NUM_THREADS"] = 12
    population_sizes = Dict(
        "influence" => Int(influence_size),
        "mutation" => Int(mutation_size),
        "cross" => Int(cross_size),
        "pop_size" => Int(mutation_size + influence_size + cross_size)
    )
    if end_rule == "generations"
        max_generations = end_number
        no_improvement_limit = Inf
    else
        if end_rule == "no_improvement"
            max_generations = Inf
            no_improvement_limit = end_number
        else
            error("Not know end rule")
        end
    end

    center_count_matrix = zeros(Int, size(M))

    non_fit_population = init_population(population_sizes["pop_size"])

    population, center_count_matrix = fitness_population(non_fit_population, center_count_matrix)
    belief_network = init_belief_network()
    belief_network = acceptance(belief_network, population, max_belief_space_size)
    best_individual = get_best(population)
    best_individuals_per_generation = Float64[]
    best_overall_individuals = Float64[]
    i = 0
    best_generation = 0
    not_improvement_count = 0
    while i < max_generations
        println("------- Current Generation ", i + 1, " -------")
        ti_cross, ti_influence, ti_mut = selection(population, population_sizes)
        ti_cross = crossover(ti_cross, crossover_type)
        ti_influence = influence(ti_influence, belief_network)
        ti_mut = mutation(ti_mut, center_count_matrix)
        population, center_count_matrix = fitness_population(vcat(ti_cross, ti_mut, ti_influence), center_count_matrix)
        belief_network = acceptance(belief_network, population, max_belief_space_size)
        best_individual_i = get_best(population)
        push!(best_individuals_per_generation, best_individual_i["obj"])

        if best_individual["obj"] > best_individual_i["obj"]
            println("New Best individual in gen ", i + 1, ". Score: ", best_individual["obj"])
            best_individual = best_individual_i
            best_generation = i
            not_improvement_count = 0
        else
            not_improvement_count += 1
        end
        push!(best_overall_individuals, best_individual["obj"])

        if not_improvement_count == no_improvement_limit
            break
        end
        

        i += 1
    end
    #identification="$(best_individual["obj"])_$(best_generation)_$(cross_size)_$(influence_size)_$(mutation)_$(size)_$(end_number)_$(max_belief_space_size)_$(crossover_type)_$(end_rule).png"
    #create_and_save_plot(best_overall_individuals, "overall_best_$(identification)")    
    #create_and_save_plot(best_individuals_per_generation, "generation_best_$(identification)")    
    return best_individual, best_generation + 1, best_individuals_per_generation, best_overall_individuals
end
