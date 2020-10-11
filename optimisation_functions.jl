## Function Definitions ##
## System Generation
function generate_system(setup)
    "Summary:   Generates the system required for the descent algorithm.

    Requires:   setup (in BabySFC/Generated_Fields, PSI/Generated_Fields,
                BabySFC/COMSOL_Fields, PSI/COMSOL_Fields).
                num_fields, m, mg, Bgoals
    Returns:    Field generating system. For the generated fields this is a list of objects
                (selection_sl (vector of vectors), selection_slc (vector of vectors),
                Bgoals (vector of functions), g (lightgraphs object), vertex_positions (vector of vectors))
                for the COMSOL fields this is a list of objects
                (output_full (LinearNDInterpolator), output_point_9 (LinearNDInterpolator))
    Code Sample Usage:
    setup = BabySFC/Generated_Fields

    num_fields = 8
    hom_coeff = 50 #coupling coefficient for the homogeneous fields P1-P3 [muT]
    grad_coeff = 50 #coupling coefficient for the linear fields P4-P8 [muT]
    m = hom_coeff*10^-6 # teslas: m rescaled for [T]
    mg = grad_coeff*10^-6 # teslas: mg rescaled for [T]

    # the goal field being targeted by the SFC
    Bgoals = [
        x -> [m, 0, 0],
        x -> [0, m, 0],
        x -> [0, 0, m],
        x -> [mg * x[1], 0,         -mg * x[3]],
        x -> [mg * x[2], mg * x[1],         0],
        x -> [0,         mg * x[2], -mg * x[3]],
        x -> [mg * x[3], 0,          mg * x[1]],
        x -> [0,         mg * x[3],  mg * x[2]]
    ]
    simpleloops, simpleloopscurrents, Bgoals, g, vertex_positions = generate_system(setup)"

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


function coupling_poi(poi)
    "Summary:   Coupling matrix of the field to the coils at position poi,
                depending on the setup this can use an active calculation of the
                field at the point of refer to a linear interpolation of the calculated
                COMSOL fields
    Requires:   simpleloops, simpleloopscurrents, Bgoals, g, vertex_positions, setup, fields_full, fields_point_9
                takes argument: (poi, vector of length divisible by 3, composition: x,y,z,x,y,z...)
    Returns:    Array{Float64,2}(length(poi), length(fields_full))"
    global simpleloops, simpleloopscurrents, Bgoals, g, vertex_positions, setup, fields_full, fields_point_9
    poi_mat = transpose(reshape(poi, (3,:)))
    poi_val = Array{Float64,2}(undef, length(poi), length(Bgoals))
    if setup == "BabySFC/Generated_Fields" || setup == "PSI/Generated_Fields"
        poi_mat = transpose(reshape(poi, (3,:)))
        poi_val = Array{Float64,2}(undef, length(poi), length(Bgoals))
        @threads for i in 1:length(Bgoals)
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
        @threads for i in 1:length(fields_full)
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
    "Summary:   Loads field configurations from .csv files and linearly interpolates between
                measured points to create a continuous field.
    Requires:   setup, hom_coeff, grad_coeff, num_fields, the corresponding field files. The
                files must be stored in the following way:
                *setup*/P*i*Field_*hom_coeff*_*grad_coeff*_Full.csv
                with file[:,1:3] x,y,z coordinates, file[:,4:6] x,y,z field components.
    Returns:    output_full (LinearNDInterpolator), output_point_9 (LinearNDInterpolator)"
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
    "Summary:   Gives the search area limits for the various field configurations.
    Requires:   setup
    Returns:    limits_inner, limits_outer, centre"
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
    limits_inner = [(MSR_dim/2 .+MSR_centre) (-(MSR_dim/2 .-MSR_centre))]
    limits_outer = [(MSR_dim_outer/2 .+MSR_centre) (-(MSR_dim_outer/2 .-MSR_centre))]
    centre = MSR_centre
    return limits_inner, limits_outer, centre
end


## Define the search limits of the algorithm
function search_area(num_sens)
    "Summary:   Gives the search area for the setup
    Requires:   setup, num_sens
    Returns:    max_vec_outer (vector length 3*num_sens),
                min_vec_outer (vector length 3*num_sens),
                max_vec_inner (vector length 3*num_sens),
                min_vec_inner (vector length 3*num_sens),
                search_centre (vector length 3*num_sens)"
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
    "Summary:   Gives a random point in the search domain.
    Requires:   setup, constrained (True/False), num_sens
    Returns:    vector length 3*num_sens"
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
    "Summary:   Loss function for the descent, returns the condition number
                of the coupling matrix for a poi.
    Requires:   poi (vector length 3*num_sens)
    Returns:    Float64"
    return cond(coupling_poi(poi))
end


## Determines the step size for the k-th step

function a_k(k)
    "Summary:   Step length for the SPSA algorithm as a function of the iteration k.
                a: initial step size, A: max iterations/10, k: iteration, alpha: step learning rate,
                c: initial perturbation size(used to scale the approximation
                to the gradient. It needs to be >0.), gamma: gradient learning rate,
                theta: previous coordinate, n: max iterations
    Requires:   a(Float), A(Float), alpha(Float), k(Int)
    Returns:    gamma(Float)"
    global a, A, alpha
    gamma = a/((k+A)^alpha)
    return gamma
end


## Sets the values to be within the search area
function constrainer(theta)
    "Summary:   Restricts the poi determined by the algorithm to be within the search domain.
    Requires:   constrained (true/false), theta(vector length num_sens*3)
    Returns:    theta(vector length num_sens*3)"
    global constrained
    if constrained == true
        max_vec_outer, min_vec_outer, max_vec_inner, min_vec_inner, search_centre = search_area(Int64(length(theta)/3))
        theta_const = max.(min_vec_outer[1:length(theta)], (min.(max_vec_outer[1:length(theta)], theta)))
        theta_const = transpose(reshape(theta_const, (3,:)))
        for i in 1:Int64(length(theta_const)/3)
            if (min_vec_inner[1] <theta_const[i,1] < max_vec_inner[1]) & (min_vec_inner[2] <theta_const[i,2] < max_vec_inner[2]) & (min_vec_inner[3] <theta_const[i,3] < max_vec_inner[3])
                for j in 1:3
                    if (theta_const[i,j] > search_centre[j])
                        theta_const[i,j] = max_vec_inner[j]
                    else
                        theta_const[i,j] = min_vec_inner[j]
                    end
                end
            end
        end
        return vec(transpose(theta_const))
    else
        return theta
    end
end

## Calculates the gradient at point theta with an SPSA method
function SPSA_gradient(theta, k)
    "Summary:   Calculates the gradient at point theta for the iteration k. The gradient
                is calculated in a random direction (given by a random vector drawn from
                a length(theta) dimensional binomial distribution). The random perturbation
                vector is set to decrease in magnitude for increasing k.
                a: initial step size, A: max iterations/10, k: iteration, alpha: step learning rate,
                c: initial perturbation size(used to scale the approximation
                to the gradient. It needs to be >0.), gamma: gradient learning rate,
                theta: previous coordinate, n: max iterations
    Requires:   theta(vector length num_sens*3), k (Int), c (Float), gamma (Float)
                d (Probability Object (Binomial(1,0.5)))
    Returns:    g_hat (vector length num_sens*3)"
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
    "Summary:   Creates the first two steps of the SPSA algorithm.
    Requires:   theta_0, theta_1 vectors of length num_sens*3
    Returns:    theta_list (vector of vectors), g_hat_list (vector of vectors)"
    theta_list = [constrainer(theta_0)]
    push!(theta_list, constrainer(theta_1))
    g_hat_list = [SPSA_gradient(theta_list[1], 1)]
    push!(g_hat_list, SPSA_gradient(theta_list[2], 2))
    return theta_list, g_hat_list
end


## Implementation of the SPSA Descent Algorithm
function SPSA(theta_0, theta_1)
    "Summary:   Conducts the SPSA gradient descent for 2 starting positions.
                a: initial step size, A: max iterations/10, k: iteration, alpha: step learning rate,
                c: initial perturbation size(used to scale the approximation
                to the gradient. It needs to be >0.), gamma: gradient learning rate,
                theta: previous coordinate, n: max iterations
    Requires:   theta_0, theta_1 (vectors of length num_sens*3)
    Returns:    theta_list (vector of vectors), cond_list (vector of vectors)"

    theta_list, g_hat_list = SPSA_initialise(theta_0, theta_1)
    cond_list = []
    for k in 2:n
        g_hat = SPSA_gradient(theta_list[k],k)
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
    "Summary:   Conducts the SPSA descent for a given number of sensors.
    Requires:   poi_0, poi_1 (vectors of length num_sens*3), num_sensors (Int)
    Returns:    poi_new (vector of vectors), conds_new (vector of vectors)"
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
    "Summary:   Turns the vector of vectors object into a 2D matrix.
    Requires:   vec_of_vec (Vectors of vectors: each length(vector) = (3*num_sens)
    Returns:    output Array{Float64,2}(length(vec_of_vec)*len,3)"
    len = Int64(length(vec_of_vec[1])/3)
    output = Array{Float64,2}(undef, length(vec_of_vec)*len,3)
    for i in 1:length(vec_of_vec)
        output[len*(i-1)+1:len*i,:] = transpose(reshape(vec_of_vec[i], (3,:)))
    end
    return output
end


function saveresult(filename, file_poi, file_cond)
    "Summary:   Saves the position (vector of vectors) and condition list to the location
                Output(setup)/poi/(filename)_(filename_suffix).csv
    Requires:   setup, filename_suffix, filename (String)
                file_poi, file_cond (Vectors of vectors: each length(vector) = (3*num_sens)
    Returns:    None"
    global setup, filename_suffix
    writedlm("Output$(setup)/poi/$(filename)_$(filename_suffix).csv", make_arrays(file_poi),  ',')
    writedlm("Output$(setup)/cond/$(filename)_$(filename_suffix).csv", file_cond,  ',')
end

# function to calculate corners given inner and outer search limits
function perm(n, mins, maxes)
    "Summary:   Gives the corners of a cube with side length
                midway between mins and maxes.
    Requires:   n (Int), mins (vector length 3), maxes (vector length 3)
    Returns:    vector length 24"
    return vec(collect(Iterators.product(collect(zip(mins,maxes))...)))
end
