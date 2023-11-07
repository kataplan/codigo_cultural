using JuMP, Gurobi

#Contador infactibles
global cont_infactible = 0;
global GRB_ENV = Gurobi.Env()

function Gurobi_optimal(C)
    num_stations = length(ESTACIONES)
    num_candidatas = length(CANDIDATAS)

    E = zeros(Int64, num_stations)
    m = Model(Gurobi.Optimizer) #No muestra resultados por consola.
    set_optimizer_attribute(m, "OutputFlag", 0)
    @variable(m, x[i=1:num_stations, j=1:num_candidatas], Bin)

    if instancia == 0
        @objective(m, Min, sum(dist[i, j] * x[i, j] for i in ESTACIONES, j in CANDIDATAS))
        for i in ESTACIONES
            @constraint(m, sum(x[i, j] for j in CANDIDATAS) == 1)
        end

        for i in ESTACIONES
            for j in CANDIDATAS
                @constraint(m, x[i, j] <= c[i, j] * C[j])
            end
        end

        for j in CANDIDATAS
            if C[j] == 1
                if allprior == 1
                    for l in PRIORIDADES
                        pxsum = @expression(m, sum(prior[i, l] * x[i, j] for i in ESTACIONES))
                        psum = @expression(m, sum(prior[i, l] for i in ESTACIONES))
                        @constraint(m, (pxsum - floor(psum / cl) * C[j]) <= prioridad)
                        @constraint(m, (floor(psum / cl) * C[j] - pxsum) <= prioridad)
                    end
                else
                    l = 1
                    pxsum = @expression(m, sum(prior[i, l] * x[i, j] for i in ESTACIONES))
                    psum = @expression(m, sum(prior[i, l] for i in ESTACIONES))
                    @constraint(m, (pxsum - floor(psum / cl) * C[j]) <= prioridad)
                    @constraint(m, (floor(psum / cl) * C[j] - pxsum) <= prioridad)
                end
            end
        end

        for j in CANDIDATAS
            expr2 = @expression(m, sum(r_menos[i] * x[i, j] for i in ESTACIONES))
            expr3 = @expression(m, sum(r_mas[i] * x[i, j] for i in ESTACIONES))
            @constraint(m, (expr2 - expr3) <= balance * (expr2 + expr3))
            @constraint(m, (-expr2 + expr3) <= balance * (expr2 + expr3))
        end


        optimize!(m)
        status = termination_status(m)
        if status == MOI.INFEASIBLE_OR_UNBOUNDED
            println("INFACTIBLE, BUSCANDO OTRA SOLUCIÓN")
            #cont_infactible = cont_infactible + 1;
            return Inf, zeros(num_stations, num_stations), [1000000]
        end

        #println("TERMINATION STATUS: ", status);
        Z_opt = objective_value(m)
        #println("FUNCION OBJETIVO POR RETORNAR: ", round(Z_opt));
        x_opt = value.(x)
        #println("X_OPT: ", length(x_opt));

        if (status != MOI.OPTIMAL && status != MOI.LOCALLY_SOLVED) || (round(Z_opt) - floor(round(Z_opt)) != 0 || length(x_opt) == 0)
            return Inf, zeros(num_stations, num_stations), [1000000]
        else
            # Create an array to store the slack values
            slack_values = zeros(length(CANDIDATAS))

            # Calculate the slack values
            for j = 1:length(CANDIDATAS)
                expr2 = @expression(m, sum(r_menos[i] * x[i, j] for i = 1:length(ESTACIONES)))
                expr3 = @expression(m, sum(r_mas[i] * x[i, j] for i = 1:length(ESTACIONES)))
                slack_values[j] = JuMP.value(expr2 - expr3) - balance * (JuMP.value(expr2) + JuMP.value(expr3))
            end

            # Create an array to store matching values per row
            matching_values = []

            # Iterate through each row of M
            for row in eachrow(M)
                slack_value = 0
                # Get the indices of centers in the same cluster                
                for i in 1:length(row)
                    if row[i] == 1 && slack_values[i] != 0
                        slack_value = slack_values[i]
                        break
                    end
                end

                # Add the values of the current row to the main array
                push!(matching_values, slack_value)
            end


            for i in ESTACIONES
                arrayE = findall(x -> x == 1, x_opt[i, :])
                if length(arrayE) > 0
                    E[i] = arrayE[1]
                end
            end
            return round(Z_opt), E, matching_values
        end
    else
        @objective(m, Min, sum(dist[i, j] * x[i, j] for i = 1:length(ESTACIONES), j = 1:length(CANDIDATAS)))
        for i = 1:length(ESTACIONES)
            @constraint(m, sum(x[i, j] for j = 1:length(CANDIDATAS)) == 1)
        end

        for i = 1:length(ESTACIONES)
            for j = 1:length(CANDIDATAS)
                @constraint(m, x[i, j] <= c[i, j] * C[j])
            end
        end

        for j = 1:length(CANDIDATAS)
            if C[j] == 1
                if allprior == 1
                    for l in PRIORIDADES
                        pxsum = @expression(m, sum(prior[i, l] * x[i, j] for i = 1:length(ESTACIONES)))
                        psum = @expression(m, sum(prior[i, l] for i = 1:length(ESTACIONES)))
                        @constraint(m, (pxsum - floor(psum / cl) * C[j]) <= prioridad)
                        @constraint(m, (floor(psum / cl) * C[j] - pxsum) <= prioridad)
                    end
                else
                    l = 1
                    pxsum = @expression(m, sum(prior[i, l] * x[i, j] for i = 1:length(ESTACIONES)))
                    psum = @expression(m, sum(prior[i, l] for i = 1:length(ESTACIONES)))
                    @constraint(m, (pxsum - floor(psum / cl) * C[j]) <= prioridad)
                    @constraint(m, (floor(psum / cl) * C[j] - pxsum) <= prioridad)
                end
            end
        end


        for j = 1:length(CANDIDATAS)
            expr2 = @expression(m, sum(r_menos[i] * x[i, j] for i = 1:length(ESTACIONES)))
            expr3 = @expression(m, sum(r_mas[i] * x[i, j] for i = 1:length(ESTACIONES)))
            @constraint(m, (expr2 - expr3) <= balance * (expr2 + expr3))
            @constraint(m, (-expr2 + expr3) <= balance * (expr2 + expr3))

        end

        optimize!(m)
        # Create an array to store the slack values



        status = termination_status(m)
        if status == MOI.INFEASIBLE_OR_UNBOUNDED
            #ultimo solución objetivo
            last_obje = 0
            println("INFACTIBLE, BUSCANDO OTRA SOLUCIÓN")
            if cont_infactible == 0
                global cont_infactible = cont_infactible + 1
                return Inf, zeros(num_stations, num_stations), [1000000]
            elseif (cont_infactible == 0 || cont_infactible > 0) && last_obje != 0
                global cont_infactible = cont_infactible + 1
                return last_obje, zeros(num_stations, num_stations), [1000000]
            else
                global cont_infactible = cont_infactible + 1
                return Inf, zeros(num_stations, num_stations), [1000000]
            end
        end
        slack_values = zeros(length(CANDIDATAS))

        # Calculate the slack values
        for j = 1:length(CANDIDATAS)
            expr2 = @expression(m, sum(r_menos[i] * x[i, j] for i = 1:length(ESTACIONES)))
            expr3 = @expression(m, sum(r_mas[i] * x[i, j] for i = 1:length(ESTACIONES)))
            slack_values[j] = JuMP.value(expr2 - expr3) - balance * (JuMP.value(expr2) + JuMP.value(expr3))
        end

        # Create an array to store matching values per row
        matching_values = []

        # Iterate through each row of M
        for row in eachrow(M)
            slack_value = 0
            # Get the indices of centers in the same cluster                
            for i in 1:length(row)
                if row[i] == 1 && slack_values[i] != 0
                    slack_value = slack_values[i]
                    break
                end
            end

            # Add the values of the current row to the main array
            push!(matching_values, slack_value)
        end
        #println("TERMINATION STATUS: ", status);
        Z_opt = objective_value(m)
        #println("FUNCION OBJETIVO POR RETORNAR: ", round(Z_opt));
        x_opt = value.(x)
        global last_obj = round(Z_opt)
        last_obje = last_obj
        #println("X_OPT: ", length(x_opt));
        # Create an array to store the slack values
        slack_values = zeros(length(CANDIDATAS))

        # Calculate the slack values
        for j = 1:length(CANDIDATAS)
            expr2 = @expression(m, sum(r_menos[i] * x[i, j] for i = 1:length(ESTACIONES)))
            expr3 = @expression(m, sum(r_mas[i] * x[i, j] for i = 1:length(ESTACIONES)))
            slack_values[j] = JuMP.value(expr2 - expr3) - balance * (JuMP.value(expr2) + JuMP.value(expr3))
        end

        # Create an array to store matching values per row
        matching_values = []

        # Iterate through each row of M
        for row in eachrow(M)
            slack_value = 0
            # Get the indices of centers in the same cluster                
            for i in 1:length(row)
                if row[i] == 1 && slack_values[i] != 0
                    slack_value = slack_values[i]
                    break
                end
            end

            # Add the values of the current row to the main array
            push!(matching_values, slack_value)
        end


        # Luego puedes usar estos valores en tu algoritmo cultural
        if (status != MOI.OPTIMAL && status != MOI.LOCALLY_SOLVED) || (round(Z_opt) - floor(round(Z_opt)) != 0 || length(x_opt) == 0)
            return Inf, zeros(num_stations, num_stations), [10000000]
        else
            slack_values = zeros(length(CANDIDATAS))

            # Calculate the slack values
            for j = 1:length(CANDIDATAS)
                expr2 = @expression(m, sum(r_menos[i] * x[i, j] for i = 1:length(ESTACIONES)))
                expr3 = @expression(m, sum(r_mas[i] * x[i, j] for i = 1:length(ESTACIONES)))
                slack_values[j] = JuMP.value(expr2 - expr3) - balance * (JuMP.value(expr2) + JuMP.value(expr3))
            end

            # Create an array to store matching values per row
            matching_values = []

            # Iterate through each row of M
            for row in eachrow(M)
                slack_value = 0
                # Get the indices of centers in the same cluster                
                for i in 1:length(row)
                    if row[i] == 1 && slack_values[i] != 0
                        slack_value = slack_values[i]
                        break
                    end
                end

                # Add the values of the current row to the main array
                push!(matching_values, slack_value)
            end
            for i = 1:length(ESTACIONES)
                arrayE = findall(x -> x == 1, x_opt[i, :])
                if length(arrayE) > 0
                    E[i] = arrayE[1]
                end
            end
            return round(Z_opt), E, matching_values
        end
    end
end
