using DelimitedFiles, LinearAlgebra, Random, Distributions, PyPlot, JLD
using PyCall
@pyimport scipy.interpolate as si
using Coils, Coils.CoilsPlot

include("optim_functions_brute_scaled.jl")


## System Initialisation

setup = "BabySFC/Generated_Fields"  #Optimisation setup: can be "BabySFC/Generated_Fields" or
                                    #"BabySFC/COMSOL_Fields" or "PSI/Generated_Fields" or "PSI/COMSOL_Fields"
constrained = true
num_fields = 8
num_sens = 8
max_vec_outer, min_vec_outer, max_vec_inner, min_vec_inner, search_centre = search_area(num_sens)
prob = 0.5
d = Binomial(1,prob)
m = 50e-6 # teslas
mg = 20e-6 # teslas

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
if setup == "BabySFC/Generated_Fields" || setup == "PSI/Generated_Fields"
    simpleloops, simpleloopscurrents, Bgoals, g, vertex_positions = generate_system(setup)
elseif setup == "BabySFC/COMSOL_Fields" || setup == "PSI/COMSOL_Fields"
    fields_full, fields_point_9 = generate_system(setup)
end

print("Determining the ideal number of sensors, ($(setup)):")
## Descent Parameters
n, a, c, alpha, gamma_var = 200, 1e-3, 0.01, 0.602, 0.101
A = n/10

initial_0 = initial_poi(3)
initial_1 = initial_poi(3)
coupling_poi(initial_0)
for num_sensors in 4:20
    print(num_sensors, ": ")
    global initial_0, initial_1
    intermediate_poi, intermediate_grads, intermediate_conds = next_num_sens(initial_0, initial_1, num_sensors)
    initial_0 = intermediate_poi[end-1]
    initial_1 = intermediate_poi[end]
    writedlm("Output$(setup)/num_sens_vals/poi_num_sens_$(num_sensors)_scaled.csv", make_arrays(intermediate_poi),  ',')
    writedlm("Output$(setup)/num_sens_vals/grads_num_sens_$(num_sensors)_scaled.csv",make_arrays(intermediate_grads),  ',')
    writedlm("Output$(setup)/num_sens_vals/conds_num_sens_$(num_sensors)_scaled.csv", intermediate_conds,  ',')
end


## New variant of the num_sensors calculation, no use of old positions
n, a, c, alpha, gamma_var = 200, 1e-3, 0.01, 0.602, 0.101
A = n/10

using Distributed
for num_sensors in 4:20
    print(num_sensors, ": ")
    intermediate_conds = [1000]
    intermediate_poi = []
    for i in 1:10
        intermediate_poi_i, intermediate_grads_i, intermediate_conds_i = next_num_sens(initial_poi(num_sensors), initial_poi(num_sensors), num_sensors)
        if last(intermediate_conds_i) < last(intermediate_conds)
            intermediate_conds = intermediate_conds_i
            intermediate_poi = intermediate_poi_i
        end
    end
    writedlm("Output$(setup)/num_sens_vals/poi_num_sens_$(num_sensors)_resetting_2.csv", make_arrays(intermediate_poi),  ',')
    writedlm("Output$(setup)/num_sens_vals/conds_num_sens_$(num_sensors)_resetting_2.csv", intermediate_conds,  ',')
end



print("Randomly Initialised Starting Conditions, Constrained, ($(setup)) ")
num_sensors = 8
n = 100
A = n/10
a = 1e-3
c = 0.01
alpha = 0.602
gamma_var = 0.101
airss_poi, airss_grad, airss_cond = [], [], []
min_poi = initial_poi(num_sensors)
min_cond = f(min_poi)
for i in 1:100 #try @parallel
    global min_cond, min_poi, min_grad
    poi_0, poi_1 = initial_poi(num_sensors), initial_poi(num_sensors)
    airss_poi_i, airss_grad_i, airss_cond_i = SPSA(poi_0, poi_1)
    push!(airss_poi, (last(airss_poi_i), airss_poi_i[1]))
    push!(airss_cond, (last(airss_cond_i), airss_cond_i[1]))
    print("Stepped: ",i, " Condition Number: ", last(airss_cond_i)," ")
    if last(min_cond) > last(airss_cond_i)
        min_cond, min_poi = airss_cond_i, airss_poi_i
        print("The minimal condition number has decreased to: " , last(min_cond), " ")
    end
end

writedlm("Output$(setup)/airss_end_conds_scaled.csv", airss_cond,  ',')

## Fine tuning of the AIRSS method
n, a, c, alpha, gamma_var = 1000, 1e-3, 0.01, 0.602, 0.101
A = n/10

poi_fine_tuned_SPSA, grad_fine_tuned_SPSA, cond_fine_tuned_SPSA = SPSA((min_poi[end-1]),(min_poi[end]))

writedlm("Output$(setup)/poi_fine_tuned_SPSA_scaled.csv", make_arrays(poi_fine_tuned_SPSA),  ',')

writedlm("Output$(setup)/cond_fine_tuned_SPSA_scaled.csv",cond_fine_tuned_SPSA,  ',')

print("Descent Complete($(setup))")


## tims poi

tims_poi = [0.5, -0.47, 0.5, 0.5, -0.47, -0.5, -0.5, -0.47, 0.5, -0.5, -0.47, -0.5, 0,1,0.5,0,1,-0.5, 0,0.315,0.485,0,0.315,-0.485]

function perm(n, mins, maxes)
    return [Tuple([(mins[i], maxes[i])[(t >> i) % 2+1] for i in 1:n]) for t in 1:(2^n)]
end

max_mat = transpose(reshape(max_vec_outer[1:3], (3,:)))
min_mat = transpose(reshape(min_vec_outer[1:3], (3,:)))
max_mat2 = transpose(reshape(max_vec_inner[1:3], (3,:)))
min_mat2 = transpose(reshape(min_vec_inner[1:3], (3,:)))


side_poi_0 = collect(Iterators.flatten(perm(3, min_mat, max_mat)))
side_poi_1 =collect(Iterators.flatten(perm(3, (min_mat+min_mat2)/2, (max_mat+max_mat2)/2)))
n, a, c, alpha, gamma_var = 1000, 1e-3, 0.01, 0.602, 0.101
A = n/10
poi_tim_tuned_SPSA, grad_tim_tuned_SPSA, cond_tim_tuned_SPSA = SPSA(initial_poi(8),tims_poi)
poi_corner_tuned_SPSA, grad_corner_tuned_SPSA, cond_corner_tuned_SPSA = SPSA(side_poi_0,side_poi_1)

writedlm("Output$(setup)/poi_tim_tuned_SPSA_COMSOL_5020_fine.csv", make_arrays(poi_tim_tuned_SPSA),  ',')

writedlm("Output$(setup)/cond_tim_tuned_SPSA_COMSOL_5020_fine.csv",cond_tim_tuned_SPSA,  ',')

writedlm("Output$(setup)/poi_corners_tuned_SPSA_COMSOL_5020_fine.csv", make_arrays(poi_corner_tuned_SPSA),  ',')

writedlm("Output$(setup)/cond_corners_tuned_SPSA_COMSOL_5020_fine.csv",cond_corner_tuned_SPSA,  ',')

print("Descent Complete($(setup))")

coupling_poi(initial_poi(8))
