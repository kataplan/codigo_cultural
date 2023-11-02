function update_center_counts(individual, count_matrix)
    for i in 1:size(individual, 1)
        if individual[i] == 1
            count_matrix[:, i] .+= M[:, i]
        end
    end
    return count_matrix
end

function calculate_variety_array(individual, normalized_matrix)
    array = zeros(Float64, size(M,1))
    mask = binary_to_mask(individual)
    for (index, value) in enumerate(mask)
        array[index] = normalized_matrix[index, value]
    end
    return array
end

function normalize_counts(count_matrix)
    normalized_matrix = zeros(Float64, size(count_matrix))

    for i in 1:size(count_matrix, 1)
        max_value = maximum(count_matrix[i, :])
        if max_value > 0
            normalized_matrix[i, :] .= count_matrix[i, :] ./ max_value
        end
    end

    return normalized_matrix
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

function select_random_position(arr, value::Int)
    valid_positions = findall(x -> x == 1 && x != value, arr)
    if isempty(valid_positions)
        return value  # No se encontraron posiciones válidas
    end

    random_position = valid_positions[rand(1:length(valid_positions))]
    return random_position
end

function check_centers(arr)
    println(count(isequal(1), arr))
    return count(isequal(1), arr) == 15
end