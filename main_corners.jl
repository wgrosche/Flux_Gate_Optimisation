using DelimitedFiles, LinearAlgebra, Random, Distributions, PyPlot, JLD
using PyCall
@pyimport scipy.interpolate as si
using Coils, Coils.CoilsPlot


include("optimisation_functions.jl")


## System Initialisation

setup = "BabySFC/Generated_Fields"  #Optimisation setup: can be "BabySFC/Generated_Fields" or
                                    #"BabySFC/COMSOL_Fields" or "PSI/Generated_Fields" or "PSI/COMSOL_Fields"
constrained = true
num_fields = 8
num_sens = 8
max_vec_outer, min_vec_outer, max_vec_inner, min_vec_inner, search_centre = search_area(num_sens)
prob = 0.5
d = Binomial(1,prob)


## Input fields homogeneous & gradient
hom_coeff = 50
grad_coeff = 5
m = hom_coeff*10^-6 # teslas
mg = grad_coeff*10^-6 # teslas
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

# generation of system based on the setup variable
if setup == "BabySFC/Generated_Fields" || setup == "PSI/Generated_Fields"
    simpleloops, simpleloopscurrents, Bgoals, g, vertex_positions = generate_system(setup)
elseif setup == "BabySFC/COMSOL_Fields" || setup == "PSI/COMSOL_Fields"
    fields_full, fields_point_9 = generate_system(setup)
end

## Savedata Options

filename_suffix = "$(hom_coeff)_$(grad_coeff)_finalised_run"
# poi calculated by tim


max_mat = transpose(reshape(max_vec_outer[1:3], (3,:)))
min_mat = transpose(reshape(min_vec_outer[1:3], (3,:)))
max_mat2 = transpose(reshape(max_vec_inner[1:3], (3,:)))
min_mat2 = transpose(reshape(min_vec_inner[1:3], (3,:)))

side_poi_0 = collect(Iterators.flatten(perm(3, min_mat, max_mat)))
side_poi_1 =collect(Iterators.flatten(perm(3, (min_mat+min_mat2)/2, (max_mat+max_mat2)/2)))

n, a, c, alpha, gamma_var = 1000, 1e-3, 0.01, 0.602, 0.101
A = n/10


corner_poi, corner_cond = SPSA(side_poi_0,side_poi_1)
saveresult("Corners", corner_poi, corner_cond )

print("Descent Complete($(setup))")



setup = "BabySFC/Generated_Fields"  #Optimisation setup: can be "BabySFC/Generated_Fields" or
                                    #"BabySFC/COMSOL_Fields" or "PSI/Generated_Fields" or "PSI/COMSOL_Fields"
constrained = true
num_fields = 8
num_sens = 8
max_vec_outer, min_vec_outer, max_vec_inner, min_vec_inner, search_centre = search_area(num_sens)
prob = 0.5
d = Binomial(1,prob)


## Input fields homogeneous & gradient
hom_coeff = 50
grad_coeff = 50
m = hom_coeff*10^-6 # teslas
mg = grad_coeff*10^-6 # teslas
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

# generation of system based on the setup variable
if setup == "BabySFC/Generated_Fields" || setup == "PSI/Generated_Fields"
    simpleloops, simpleloopscurrents, Bgoals, g, vertex_positions = generate_system(setup)
elseif setup == "BabySFC/COMSOL_Fields" || setup == "PSI/COMSOL_Fields"
    fields_full, fields_point_9 = generate_system(setup)
end

## Savedata Options

filename_suffix = "$(hom_coeff)_$(grad_coeff)_finalised_run"
# poi calculated by tim


max_mat = transpose(reshape(max_vec_outer[1:3], (3,:)))
min_mat = transpose(reshape(min_vec_outer[1:3], (3,:)))
max_mat2 = transpose(reshape(max_vec_inner[1:3], (3,:)))
min_mat2 = transpose(reshape(min_vec_inner[1:3], (3,:)))

side_poi_0 = collect(Iterators.flatten(perm(3, min_mat, max_mat)))
side_poi_1 =collect(Iterators.flatten(perm(3, (min_mat+min_mat2)/2, (max_mat+max_mat2)/2)))

n, a, c, alpha, gamma_var = 1000, 1e-3, 0.01, 0.602, 0.101
A = n/10


corner_poi, corner_cond = SPSA(side_poi_0,side_poi_1)
saveresult("Corners", corner_poi, corner_cond )

print("Descent Complete($(setup))")



setup = "BabySFC/COMSOL_Fields"  #Optimisation setup: can be "BabySFC/Generated_Fields" or
                                    #"BabySFC/COMSOL_Fields" or "PSI/Generated_Fields" or "PSI/COMSOL_Fields"
constrained = true
num_fields = 8
num_sens = 8
max_vec_outer, min_vec_outer, max_vec_inner, min_vec_inner, search_centre = search_area(num_sens)
prob = 0.5
d = Binomial(1,prob)


## Input fields homogeneous & gradient
hom_coeff = 50
grad_coeff = 5
m = hom_coeff*10^-6 # teslas
mg = grad_coeff*10^-6 # teslas
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

# generation of system based on the setup variable
if setup == "BabySFC/Generated_Fields" || setup == "PSI/Generated_Fields"
    simpleloops, simpleloopscurrents, Bgoals, g, vertex_positions = generate_system(setup)
elseif setup == "BabySFC/COMSOL_Fields" || setup == "PSI/COMSOL_Fields"
    fields_full, fields_point_9 = generate_system(setup)
end

## Savedata Options

filename_suffix = "$(hom_coeff)_$(grad_coeff)_finalised_run"
# poi calculated by tim


max_mat = transpose(reshape(max_vec_outer[1:3], (3,:)))
min_mat = transpose(reshape(min_vec_outer[1:3], (3,:)))
max_mat2 = transpose(reshape(max_vec_inner[1:3], (3,:)))
min_mat2 = transpose(reshape(min_vec_inner[1:3], (3,:)))

side_poi_0 = collect(Iterators.flatten(perm(3, min_mat, max_mat)))
side_poi_1 =collect(Iterators.flatten(perm(3, (min_mat+min_mat2)/2, (max_mat+max_mat2)/2)))

n, a, c, alpha, gamma_var = 10000, 1e-3, 0.01, 0.602, 0.101
A = n/10


corner_poi, corner_cond = SPSA(side_poi_0,side_poi_1)
saveresult("Corners", corner_poi, corner_cond )

print("Descent Complete($(setup))")


setup = "BabySFC/COMSOL_Fields"  #Optimisation setup: can be "BabySFC/Generated_Fields" or
                                    #"BabySFC/COMSOL_Fields" or "PSI/Generated_Fields" or "PSI/COMSOL_Fields"
constrained = true
num_fields = 8
num_sens = 8
max_vec_outer, min_vec_outer, max_vec_inner, min_vec_inner, search_centre = search_area(num_sens)
prob = 0.5
d = Binomial(1,prob)


## Input fields homogeneous & gradient
hom_coeff = 50
grad_coeff = 50
m = hom_coeff*10^-6 # teslas
mg = grad_coeff*10^-6 # teslas
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

# generation of system based on the setup variable
if setup == "BabySFC/Generated_Fields" || setup == "PSI/Generated_Fields"
    simpleloops, simpleloopscurrents, Bgoals, g, vertex_positions = generate_system(setup)
elseif setup == "BabySFC/COMSOL_Fields" || setup == "PSI/COMSOL_Fields"
    fields_full, fields_point_9 = generate_system(setup)
end

## Savedata Options

filename_suffix = "$(hom_coeff)_$(grad_coeff)_finalised_run"
# poi calculated by tim


max_mat = transpose(reshape(max_vec_outer[1:3], (3,:)))
min_mat = transpose(reshape(min_vec_outer[1:3], (3,:)))
max_mat2 = transpose(reshape(max_vec_inner[1:3], (3,:)))
min_mat2 = transpose(reshape(min_vec_inner[1:3], (3,:)))

side_poi_0 = collect(Iterators.flatten(perm(3, min_mat, max_mat)))
side_poi_1 =collect(Iterators.flatten(perm(3, (min_mat+min_mat2)/2, (max_mat+max_mat2)/2)))

n, a, c, alpha, gamma_var = 10000, 1e-3, 0.01, 0.602, 0.101
A = n/10


corner_poi, corner_cond = SPSA(side_poi_0,side_poi_1)
saveresult("Corners", corner_poi, corner_cond )

print("Descent Complete($(setup))")



setup = "PSI/COMSOL_Fields"  #Optimisation setup: can be "BabySFC/Generated_Fields" or
                                    #"BabySFC/COMSOL_Fields" or "PSI/Generated_Fields" or "PSI/COMSOL_Fields"
constrained = true
num_fields = 8
num_sens = 8
max_vec_outer, min_vec_outer, max_vec_inner, min_vec_inner, search_centre = search_area(num_sens)
prob = 0.5
d = Binomial(1,prob)


## Input fields homogeneous & gradient
hom_coeff = 50
grad_coeff = 5
m = hom_coeff*10^-6 # teslas
mg = grad_coeff*10^-6 # teslas
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

# generation of system based on the setup variable
if setup == "BabySFC/Generated_Fields" || setup == "PSI/Generated_Fields"
    simpleloops, simpleloopscurrents, Bgoals, g, vertex_positions = generate_system(setup)
elseif setup == "BabySFC/COMSOL_Fields" || setup == "PSI/COMSOL_Fields"
    fields_full, fields_point_9 = generate_system(setup)
end

## Savedata Options

filename_suffix = "$(hom_coeff)_$(grad_coeff)_finalised_run"
# poi calculated by tim


max_mat = transpose(reshape(max_vec_outer[1:3], (3,:)))
min_mat = transpose(reshape(min_vec_outer[1:3], (3,:)))
max_mat2 = transpose(reshape(max_vec_inner[1:3], (3,:)))
min_mat2 = transpose(reshape(min_vec_inner[1:3], (3,:)))

side_poi_0 = collect(Iterators.flatten(perm(3, min_mat, max_mat)))
side_poi_1 =collect(Iterators.flatten(perm(3, (min_mat+min_mat2)/2, (max_mat+max_mat2)/2)))

n, a, c, alpha, gamma_var = 10000, 1e-3, 0.01, 0.602, 0.101
A = n/10


corner_poi, corner_cond = SPSA(side_poi_0,side_poi_1)
saveresult("Corners", corner_poi, corner_cond )

print("Descent Complete($(setup))")



setup = "PSI/COMSOL_Fields"  #Optimisation setup: can be "BabySFC/Generated_Fields" or
                                    #"BabySFC/COMSOL_Fields" or "PSI/Generated_Fields" or "PSI/COMSOL_Fields"
constrained = true
num_fields = 8
num_sens = 8
max_vec_outer, min_vec_outer, max_vec_inner, min_vec_inner, search_centre = search_area(num_sens)
prob = 0.5
d = Binomial(1,prob)


## Input fields homogeneous & gradient
hom_coeff = 50
grad_coeff = 50
m = hom_coeff*10^-6 # teslas
mg = grad_coeff*10^-6 # teslas
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

# generation of system based on the setup variable
if setup == "BabySFC/Generated_Fields" || setup == "PSI/Generated_Fields"
    simpleloops, simpleloopscurrents, Bgoals, g, vertex_positions = generate_system(setup)
elseif setup == "BabySFC/COMSOL_Fields" || setup == "PSI/COMSOL_Fields"
    fields_full, fields_point_9 = generate_system(setup)
end

## Savedata Options

filename_suffix = "$(hom_coeff)_$(grad_coeff)_finalised_run"
# poi calculated by tim


max_mat = transpose(reshape(max_vec_outer[1:3], (3,:)))
min_mat = transpose(reshape(min_vec_outer[1:3], (3,:)))
max_mat2 = transpose(reshape(max_vec_inner[1:3], (3,:)))
min_mat2 = transpose(reshape(min_vec_inner[1:3], (3,:)))

side_poi_0 = collect(Iterators.flatten(perm(3, min_mat, max_mat)))
side_poi_1 =collect(Iterators.flatten(perm(3, (min_mat+min_mat2)/2, (max_mat+max_mat2)/2)))

n, a, c, alpha, gamma_var = 10000, 1e-3, 0.01, 0.602, 0.101
A = n/10


corner_poi, corner_cond = SPSA(side_poi_0,side_poi_1)
saveresult("Corners", corner_poi, corner_cond )

print("Descent Complete($(setup))")



setup = "PSI/Generated_Fields"  #Optimisation setup: can be "BabySFC/Generated_Fields" or
                                    #"BabySFC/COMSOL_Fields" or "PSI/Generated_Fields" or "PSI/COMSOL_Fields"
constrained = true
num_fields = 8
num_sens = 8
max_vec_outer, min_vec_outer, max_vec_inner, min_vec_inner, search_centre = search_area(num_sens)
prob = 0.5
d = Binomial(1,prob)


## Input fields homogeneous & gradient
hom_coeff = 50
grad_coeff = 5
m = hom_coeff*10^-6 # teslas
mg = grad_coeff*10^-6 # teslas
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

# generation of system based on the setup variable
if setup == "BabySFC/Generated_Fields" || setup == "PSI/Generated_Fields"
    simpleloops, simpleloopscurrents, Bgoals, g, vertex_positions = generate_system(setup)
elseif setup == "BabySFC/COMSOL_Fields" || setup == "PSI/COMSOL_Fields"
    fields_full, fields_point_9 = generate_system(setup)
end

## Savedata Options

filename_suffix = "$(hom_coeff)_$(grad_coeff)_finalised_run"
# poi calculated by tim


max_mat = transpose(reshape(max_vec_outer[1:3], (3,:)))
min_mat = transpose(reshape(min_vec_outer[1:3], (3,:)))
max_mat2 = transpose(reshape(max_vec_inner[1:3], (3,:)))
min_mat2 = transpose(reshape(min_vec_inner[1:3], (3,:)))

side_poi_0 = collect(Iterators.flatten(perm(3, min_mat, max_mat)))
side_poi_1 =collect(Iterators.flatten(perm(3, (min_mat+min_mat2)/2, (max_mat+max_mat2)/2)))

n, a, c, alpha, gamma_var = 1000, 1e-3, 0.01, 0.602, 0.101
A = n/10


corner_poi, corner_cond = SPSA(side_poi_0,side_poi_1)
saveresult("Corners", corner_poi, corner_cond )

print("Descent Complete($(setup))")




setup = "PSI/Generated_Fields"  #Optimisation setup: can be "BabySFC/Generated_Fields" or
                                    #"BabySFC/COMSOL_Fields" or "PSI/Generated_Fields" or "PSI/COMSOL_Fields"
constrained = true
num_fields = 8
num_sens = 8
max_vec_outer, min_vec_outer, max_vec_inner, min_vec_inner, search_centre = search_area(num_sens)
prob = 0.5
d = Binomial(1,prob)


## Input fields homogeneous & gradient
hom_coeff = 50
grad_coeff = 50
m = hom_coeff*10^-6 # teslas
mg = grad_coeff*10^-6 # teslas
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

# generation of system based on the setup variable
if setup == "BabySFC/Generated_Fields" || setup == "PSI/Generated_Fields"
    simpleloops, simpleloopscurrents, Bgoals, g, vertex_positions = generate_system(setup)
elseif setup == "BabySFC/COMSOL_Fields" || setup == "PSI/COMSOL_Fields"
    fields_full, fields_point_9 = generate_system(setup)
end

## Savedata Options

filename_suffix = "$(hom_coeff)_$(grad_coeff)_finalised_run"
# poi calculated by tim


max_mat = transpose(reshape(max_vec_outer[1:3], (3,:)))
min_mat = transpose(reshape(min_vec_outer[1:3], (3,:)))
max_mat2 = transpose(reshape(max_vec_inner[1:3], (3,:)))
min_mat2 = transpose(reshape(min_vec_inner[1:3], (3,:)))

side_poi_0 = collect(Iterators.flatten(perm(3, min_mat, max_mat)))
side_poi_1 =collect(Iterators.flatten(perm(3, (min_mat+min_mat2)/2, (max_mat+max_mat2)/2)))

n, a, c, alpha, gamma_var = 1000, 1e-3, 0.01, 0.602, 0.101
A = n/10


corner_poi, corner_cond = SPSA(side_poi_0,side_poi_1)
saveresult("Corners", corner_poi, corner_cond )

print("Descent Complete($(setup))")
