using DelimitedFiles, LinearAlgebra, Random, Distributions, PyPlot, JLD
using PyCall
@pyimport scipy.interpolate as si
using Coils, Coils.CoilsPlot

## Import the functions for the optimisation
include("optimisation_functions.jl")

## System Initialisation

setup = "BabySFC/Generated_Fields"  #Optimisation setup: can be "BabySFC/Generated_Fields" or
                                    #"BabySFC/COMSOL_Fields" or "PSI/Generated_Fields" or "PSI/COMSOL_Fields"

num_fields = 8 #number of fields imported (from the polynomial decomposition of the fields)
num_sens = 8 #number of flux gates being optimised

# the search domain that is being examined when constrained is given to be true
max_vec_outer, min_vec_outer, max_vec_inner, min_vec_inner, search_centre = search_area(num_sens)
constrained = true # determines whether or not the search area will be limited to an area around the mu metal shield

prob = 0.5 #used for the Bernoulli distribution below.
d = Binomial(1,prob) #Bernoulli distribution, used to generate random perturbation vectors.


## Input fields homogeneous & gradient
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

# generation of system based on the setup variable, in the case of the generated fields it returns the
# coil setup as given in the coils.jl package and the loops and corresponding currents through them.
# in the case of the COMSOL fields it loads in the fields from the field file (matrix csv) and linearly
# interpolates between the support points.
if setup == "BabySFC/Generated_Fields" || setup == "PSI/Generated_Fields"
    simpleloops, simpleloopscurrents, Bgoals, g, vertex_positions = generate_system(setup)
elseif setup == "BabySFC/COMSOL_Fields" || setup == "PSI/COMSOL_Fields"
    fields_full, fields_point_9 = generate_system(setup)
end

## Savedata Options
# gives the save file extension for the generated poi and condition lists
filename_suffix = "$(hom_coeff)_$(grad_coeff)_finalised_run"
## num_sensors calculation
print("Determining the ideal number of sensors, ($(setup)):")
# Descent Parameters: n: number of iterations for the descent, others used to determine a_k, c_k, refer to
# optimise_functions.jl for more detail.
n, a, c, alpha, gamma_var = 200, 1e-3, 0.01, 0.602, 0.101
A = n/10


# randomy initialises for each number of sensors. Each number of sensors is computed 10 times and the lowest
# condition number achieved is used.
for num_sensors in 4:20
    print(num_sensors, ": ")
    intermediate_conds = [1000]
    intermediate_poi = []
    for i in 1:10
        intermediate_poi_i, intermediate_conds_i = next_num_sens(initial_poi(num_sensors), initial_poi(num_sensors), num_sensors)
        if last(intermediate_conds_i) < last(intermediate_conds)
            intermediate_conds = intermediate_conds_i
            intermediate_poi = intermediate_poi_i
        end
    end
    writedlm("Output$(setup)/num_sens_vals/poi_$(num_sensors)_$(filename_suffix).csv", make_arrays(intermediate_poi),  ',')
    writedlm("Output$(setup)/num_sens_vals/conds_$(num_sensors)_$(filename_suffix).csv", intermediate_conds,  ',')
end


# AIRSS Descent calculation
# Each descent is randomly initialised and the lowest condition number achieved over 100 runs is
# passed to the fine tuning descent.
print("Randomly Initialised Starting Conditions, Constrained, ($(setup)) ")
# Descent Parameters
n, a, c, alpha, gamma_var = 200, 1e-3, 0.01, 0.602, 0.101
A = n/10
num_sensors = 8

airss_poi, airss_grad, airss_cond = [], [], []
min_poi = initial_poi(num_sensors)
min_cond = f(min_poi)
for i in 1:100 #try @parallel
    global min_cond, min_poi, min_grad
    poi_0, poi_1 = initial_poi(num_sensors), initial_poi(num_sensors)
    airss_poi_i, airss_cond_i = SPSA(poi_0, poi_1)
    push!(airss_poi, (last(airss_poi_i), airss_poi_i[1]))
    push!(airss_cond, (last(airss_cond_i), airss_cond_i[1]))
    print("Stepped: ",i, " Condition Number: ", last(airss_cond_i)," ")
    if last(min_cond) > last(airss_cond_i)
        min_cond, min_poi = airss_cond_i, airss_poi_i
        print("The minimal condition number has decreased to: " , last(min_cond), " ")
    end
end

writedlm("Output$(setup)/airss_end_conds_$(filename_suffix).csv", airss_cond,  ',')

# Fine tuning of the AIRSS method
n, a, c, alpha, gamma_var = 1000, 1e-3, 0.01, 0.602, 0.101
A = n/10

tuning_poi, tuning_conds = SPSA((min_poi[end-1]),(min_poi[end]))
saveresult("fine_tuned", tuning_poi, tuning_conds)


## Static Starting conditions: BabySFC Poi calculated by Tim Roethlisberger, Poi in the corners of the search domain for BabySFC and PSI

# initialised using the flux_gate positions given by Tim and fine tuned using the SPSA algorithm
# poi calculated by tim
tims_poi = [0.5, -0.47, 0.5, 0.5, -0.47, -0.5, -0.5, -0.47, 0.5, -0.5, -0.47, -0.5, 0,1,0.5,0,1,-0.5, 0,0.315,0.485,0,0.315,-0.485]

n, a, c, alpha, gamma_var = 1000, 1e-2, 0.01, 0.602, 0.101
A = n/10

tim_poi, tim_cond = SPSA(initial_poi(8),tims_poi)

saveresult("Tim", tim_poi, tim_cond )

# Starting points chosen to be the corners midway between the inner and outer
# limits of the search domain and the inner corners of the search domain.

max_mat = transpose(reshape(max_vec_outer[1:3], (3,:)))
min_mat = transpose(reshape(min_vec_outer[1:3], (3,:)))
max_mat2 = transpose(reshape(max_vec_inner[1:3], (3,:)))
min_mat2 = transpose(reshape(min_vec_inner[1:3], (3,:)))

side_poi_0 = collect(Iterators.flatten(perm(3, min_mat, max_mat)))
side_poi_1 =collect(Iterators.flatten(perm(3, (min_mat+min_mat2)/2, (max_mat+max_mat2)/2)))


n, a, c, alpha, gamma_var = 1000, 1e-2, 0.01, 0.602, 0.101
A = n/10

corner_poi, corner_cond = SPSA(side_poi_0,side_poi_1)
saveresult("Corners", corner_poi, corner_cond )

print("Descent Complete($(setup))")
