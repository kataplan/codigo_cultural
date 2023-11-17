function update_center_counts(individual, count_matrix)
    for i in 1:size(individual, 1)
        if individual[i] == 1
            count_matrix[:, i] .+= M[:, i]
        end
    end
    return count_matrix
end

function calculate_count_array(individual, variety_matrix)
    array = zeros(Float64, size(M, 1))
    for (index, value) in enumerate(individual)
        array[index] = variety_matrix[index, value]
    end
    return array
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
    new_array = copy(arr)
    new_array[value] = 0
    valid_positions = findall(x -> x == 1, new_array)
    if isempty(valid_positions)
        return value  # No se encontraron posiciones válidas
    end

    random_position = valid_positions[rand(1:length(valid_positions))]
    return random_position
end

function check_centers(arr)
    return count(isequal(1), arr) == 15
end

function create_and_save_plot(data::Vector{Float64}, filename::AbstractString)

    output_folder = "plots"

    output_filename = filename
    max_value = maximum(data)
    ylim_max = max_value * 1.1  # Ajusta el factor según tus datos
    output_path = joinpath(output_folder, output_filename)

    # Crear un gráfico de dispersión
    plot(data, seriestype=:line, legend=false, grid=false, xlabel="Generations", ylabel="Objective value")

    # Guardar el gráfico en la ubicación especificada
    savefig(output_path)

    println("Gráfico guardado en $output_path")
end