using Coils
using Coils.CoilsPlot
using DelimitedFiles
using JLD
using PyPlot
using LinearAlgebra
using Random
using PyCall
@pyimport scipy.interpolate as si
using LightGraphs

include("optimisation_functions.jl")
##investigating the stability of the condition number at the end of the optimisation

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


function pos_shift(poi,magnitude)
    cond_shift = Array{Float64,2}(undef, 8, 3)
    for flux_gate in 1:8
        for dir in 1:3
            poi_shift =  zeros(8, 3)
            poi_shift[flux_gate, dir] = magnitude
            cond_shift[flux_gate,dir] = f(poi.+poi_shift) - f(poi)
        end
    end
    return cond_shift
end

function random_shift(poi, magnitude)
    delta =  rand(d, length(poi))*2 .-1
    shift = delta*magnitude
    poi_shifted = poi + shift
    return poi_shifted
end

output = Array{Float64,2}(undef, length(range(-1,stop = 1, step = 0.001)),2)
##output
configuration = "fine_tuned" #Corners or Tim or fine_tuned
scaling = "50_5"
poi  = readdlm("Output$(setup)/poi/$(configuration)_$(scaling)_finalised_run.csv" , ',' , Float64)[end-7:end,:]


output2 = []
for i in eachindex(range(-1,stop = 1, step = 0.001))
    push!(output, f(random_shift(poi, i)))
end


for j in eachindex(range(-1,stop = 1, step = 0.001))
    output[i,:] = [(range(-1,stop = 1, step = 0.001))[i] , pos_shift(poi, (range(-1,stop = 1, step = 0.001)))]
end

writedlm("D:/UserSecondary/Polybox/Semesterarbeit//optimisation/Optimised_points/positional_stability.csv",output,  ',')

[(range(-1,stop = 1, step = 0.001))[1] findmax(pos_shift(final_poi, (range(-1,stop = 1, step = 0.001))[1]))[1]]
output
