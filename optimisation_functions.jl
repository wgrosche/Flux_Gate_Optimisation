## Function Definitions ##
## System Generation
## Baby SFC
function generate_system(setup)
    global num_fields, m, mg, Bgoals
    if setup == "BabySFC/Generated_Fields"
        tile_size = 0.262 #m
        ntiles = [5, 9, 5]
        totalsize = ntiles .* tile_size

        ## Generate POI and vertex_positions
        g, vertex_positions = cuboid_system(totalsize, ntiles,
            skipfaces = [false, false, true, false, false, false])

        #poi = cuboid_poi([1.7, 2.74, 1.7], [0.0, 0.0, 0.0], [30, 60, 30], filled = true)
        poi = cuboid_poi([1.7, 2.74, 1.7], [0.0, 0.0, 0.0], [30, 60, 30], filled = true)
        # prepare the large cells
        front_face_unordered = [ i for i in eachindex(vertex_positions)
            if vertex_positions[i][2] < -totalsize[2] / 2 + 0.01 ]
        front_face = order_cell(g, front_face_unordered)

        initialcells = [front_face]
        elemcurrents = [5.0, 1.0, 0.2]

        dirs = [1,2,3]

        ## Span for which the field is calculated
        spanA = extrema([ p[dirs[1]] for p in vertex_positions])
        spanB = extrema([ p[dirs[2]] for p in vertex_positions])
        spanC = extrema([ p[dirs[3]] for p in vertex_positions])
        surround = -0.05 #m
        range_x = spanA[2] - spanA[1] + surround
        range_y = spanB[2] - spanB[1] + surround
        range_z = spanC[2] - spanC[1] + surround
        #centre:
        centre_x = spanA[1] + (spanA[2]- spanA[1])/2
        centre_y = spanB[1] + (spanB[2]- spanB[1])/2
        centre_z = spanC[1] + (spanC[2]- spanC[1])/2
        poi2 = cuboid_poi([0.97,0.97,0.87], [0,0.315,0], [30, 60, 30], filled = true)
        simpleloops = []
        simpleloopscurrents = []

        for Bgoal in Bgoals
           simpleloopsi, simpleloopscurrentsi = solve_system(g, vertex_positions, poi2, Bgoal,
               verbose = true,
               tolerance = 0.0001,
               bigfloat = false,
               λ = 0, # regularization parameter
               mincurrent = minimum(elemcurrents) / 2,
               initialcells = initialcells)
           push!(simpleloops, simpleloopsi)
           push!(simpleloopscurrents, simpleloopscurrentsi)
        end
        return simpleloops, simpleloopscurrents, Bgoals, g, vertex_positions
    elseif setup == "PSI/Generated_Fields"
        vertex_positions = load("PSI/Generated_Fields/grid.jld", "vertex_positions")
        cells = load("PSI/Generated_Fields/grid.jld", "cells")
        g = load("PSI/Generated_Fields/grid.jld", "g")
        selection_slc = load("PSI/Generated_Fields/simplified_system.jld", "selection_slc")
        selection_sl = load("PSI/Generated_Fields/simplified_system.jld", "selection_sl")
        return selection_sl, selection_slc, Bgoals, g, vertex_positions
    elseif setup == "BabySFC/COMSOL_Fields" || setup == "PSI/COMSOL_Fields"
        return load_data(num_fields)
    end
end

"Coupling matrix of the field to the coils at position poi,
depending on the setup this can use an active calculation of the
field at the point of refer to a linear interpolation of the calculated
COMSOL fields"
function coupling_poi(poi)
    global simpleloops, simpleloopscurrents, Bgoals, g, vertex_positions, setup, fields_full, fields_point_9
    poi_mat = transpose(reshape(poi, (3,:)))
    poi_val = Array{Float64,2}(undef, length(poi), length(Bgoals))
    if setup == "BabySFC/Generated_Fields" || setup == "PSI/Generated_Fields"
        poi_mat = transpose(reshape(poi, (3,:)))
        poi_val = Array{Float64,2}(undef, length(poi), length(Bgoals))
        for i in 1:length(Bgoals)
            if i >= 4
                c = m/mg
            else
                c = 1
            end
            for j in 1:Int64(length(poi)/3)
                poi_val[3*(j-1)+1:3*j,i] = 10*c*(field_loops(poi_mat[j,:], simpleloops[i], simpleloopscurrents[i], vertex_positions)-field_loops(poi_mat[j,:], simpleloops[i], simpleloopscurrents[i]*0.9, vertex_positions))
            end
        end
    else
        poi_mat = transpose(reshape(poi, (3,:)))
        poi_val = Array{Float64,2}(undef, length(poi), length(fields_full))
        for i in 1:length(fields_full)
            if i >= 4
                c = m/mg
            else
                c = 1
            end
            poi_val[:,i] = 10*c*(vec(transpose(fields_full[i](poi_mat) - fields_point_9[i](poi_mat))))
        end
    end
    return poi_val
end




## Load Data
function load_data(num_fields)
    global setup, hom_coeff, grad_coeff
    output_full = []
    output_point_9 = []
    for i in 1:num_fields
        input_full = readdlm("$(setup)/P$(i)Field_$(hom_coeff)_$(grad_coeff)_Full.csv" , ',' , Float64)
        input_point_9 = readdlm("$(setup)/P$(i)Field_$(hom_coeff)_$(grad_coeff)_Point_9.csv" , ',' , Float64)
        P_i_field_full = si.LinearNDInterpolator(input_full[:,1:3], input_full[:,4:6])
        P_i_field_point_9 = si.LinearNDInterpolator(input_point_9[:,1:3], input_point_9[:,4:6])
        push!(output_full, P_i_field_full)
        push!(output_point_9, P_i_field_point_9)
    end
    return output_full, output_point_9
end

## Area around the Mu Metal Shield
function search_area_limits()
    global setup
    if setup == "BabySFC/Generated_Fields" || setup == "BabySFC/COMSOL_Fields"
        MSR_dim = [0.97, 0.97, 0.87] #in m
        MSR_centre = [0,0.315,0] #in m
        MSR_dim_outer = MSR_dim .+ 0.2
    elseif setup == "PSI/Generated_Fields" || setup == "PSI/COMSOL_Fields"
        MSR_dim = [5044.5, 5044.5, 4716.15]*10^-3 #in m
        MSR_centre = [0,0,143.625]*10^-3 #in m
        MSR_dim_outer = MSR_dim .+ 0.4
    end
    limits_inner = [(MSR_dim.+MSR_centre)/2 (-(MSR_dim.-MSR_centre)/2)]
    limits_outer = [(MSR_dim_outer.+MSR_centre)/2 (-(MSR_dim_outer.-MSR_centre)/2)]
    centre = MSR_centre
    return limits_inner, limits_outer, centre
end


## Define the search limits of the algorithm
function search_area(num_sens)
    global setup
    dims = num_sens*3
    limits_inner, limits_outer, centre = search_area_limits()
    max_vec_outer = Array{Float64,1}(undef, dims)
    min_vec_outer = Array{Float64,1}(undef, dims)
    max_vec_inner = Array{Float64,1}(undef, dims)
    min_vec_inner = Array{Float64,1}(undef, dims)
    search_centre = Array{Float64,1}(undef, dims)
    for i in 1:num_sens
        max_vec_outer[3(i-1)+1:3i] = limits_outer[:,1]
        min_vec_outer[3(i-1)+1:3i] = limits_outer[:,2]
        max_vec_inner[3(i-1)+1:3i] = limits_inner[:,1]
        min_vec_inner[3(i-1)+1:3i] = limits_inner[:,2]
        search_centre[3(i-1)+1:3i] = centre
    end
    return max_vec_outer, min_vec_outer, max_vec_inner, min_vec_inner, search_centre
end


## Initialisation (picks a point at random from the search domain)
function initial_poi(num_sens)
    limits_inner, limits_outer, centre = search_area_limits()
    n = 50
    Xs = range(limits_outer[1,1],stop = limits_outer[1,2], length = n)
    Ys = range(limits_outer[2,1],stop = limits_outer[2,2], length = n)
    Zs = range(limits_outer[3,1],stop = limits_outer[3,2], length = n)
    poi = Array{Float64,2}(undef, num_sens, 3)
    for i in 1:num_sens
        poi[i,:] = [rand(Xs), rand(Ys), rand(Zs)]
    end
    return constrainer(vec(poi))
end


## Chains the loss and matrix coupling functions, returns the condition number at points poi
function f(poi)
    return cond(coupling_poi(poi))
end


## Determines the step size for the k-th step

function a_k(k)
    global a, A, alpha
    gamma = a/((k+A)^alpha)
    return gamma
end


## Sets the values to be within the search area
function constrainer(theta)
    global constrained
    if constrained == true
        max_vec_outer, min_vec_outer, max_vec_inner, min_vec_inner, search_centre = search_area(length(theta))
        theta = max.(min_vec_outer[1:length(theta)], (min.(max_vec_outer[1:length(theta)], theta)))
        for i in 1:length(theta)
            if (theta[i] > min_vec_inner[i]) & (theta[i] < max_vec_inner[i])
                if (theta[i] > search_centre[i])
                    theta[i] = max_vec_inner[i]
                else
                    theta[i] = min_vec_inner[i]
                end
            end
        end
        return theta
    else
        return theta
    end
end


## Calculates the gradient at point theta with an SPSA method
function SPSA_gradient(theta, k)
    global c, gamma_var, d
    c_k = c/(k^gamma_var)
    delta =  rand(d, length(theta))*2 .-1
    theta_plus = theta + c_k*delta
    theta_minus = theta - c_k*delta
    y_plus = f(theta_plus)
    y_minus = f(theta_minus)
    g_hat = (y_plus - y_minus) ./ (2*c_k*delta)
    return g_hat
end


## Sets the starting points for the optimisation
function SPSA_initialise(theta_0, theta_1)
    theta_list = [constrainer(theta_0)]
    push!(theta_list, constrainer(theta_1))
    g_hat_list = [SPSA_gradient(theta_list[1], 1)]
    push!(g_hat_list, SPSA_gradient(theta_list[2], 2))
    return theta_list, g_hat_list
end


## Implementation of the SPSA Descent Algorithm
function SPSA(theta_0, theta_1)
    ## a: initial step size, A: max iterations/10, k: iteration, alpha: step learning rate, c: initial perturbation size(used to scale the approximation
    ##        to the gradient. It needs to be >0.), gamma: gradient learning rate, p: , theta: previous coordinate, n: max iterations
    global max_vec_outer, min_vec_outer, max_vec_inner, min_vec_inner, search_centre
    global A, c, a, alpha, gamma_var, n

    theta_list, g_hat_list = SPSA_initialise(theta_0, theta_1)
    cond_list = []
    for k in 2:n
        a_k_var = a/((k+A)^alpha)
        g_hat = SPSA_gradient(theta_list[k],k)
        #a_k_var = a_k(theta_list[k], theta_list[k-1], g_hat_list[k], g_hat_list[k-1],k)
        a_k_var = a_k(k)
        theta = theta_list[k] - a_k_var*g_hat
        theta = constrainer(theta)
        if (k % 100)==0
            print("Loss: ",f(theta),", ")
        end
        push!(cond_list, f(theta))
        push!(theta_list, theta)
        push!(g_hat_list, g_hat)
    end
    return theta_list, cond_list
end


## Used in the detemination of the ideal number of flux gates
function next_num_sens(poi_0, poi_1, num_sensors)
    p = num_sensors*3
    new_poi_1 = initial_poi(1)
    new_poi_0 = initial_poi(1)
    theta_1 = vcat(vec(poi_1), vec(new_poi_1))
    theta_0 = vcat(vec(poi_0), vec(new_poi_0))
    poi_new, conds_new = SPSA(theta_0, theta_1)
    return poi_new, conds_new
end


## Reshape the output vectors into arrays of the shape (3*num_sens,num_fields)
function make_arrays(vec_of_vec)
    len = Int64(length(vec_of_vec[1])/3)
    output = Array{Float64,2}(undef, length(vec_of_vec)*len,3)
    for i in 1:length(vec_of_vec)
        output[len*(i-1)+1:len*i,:] = reshape(vec_of_vec[i], (:,3))
    end
    return output
end


function saveresult(filename, file_poi, file_cond)
    global setup, filename_suffix
    writedlm("Output$(setup)/poi/$(filename)_$(filename_suffix).csv", make_arrays(file_poi),  ',')
    writedlm("Output$(setup)/cond/$(filename)_$(filename_suffix).csv", file_cond,  ',')
end
