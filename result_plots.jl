using PyPlot
using DelimitedFiles
using Coils
using Coils.CoilsPlot
using LightGraphs
#using PyPlot
# using Plots
# pyplot()
# Plots.PyPlotBackend()

include("plotting_functions.jl")


## Poi Plots

## BabySFC
setup = "OutputBabySFC/Generated_Fields"
file = "tims_tuned_SPSA_scaled"

titles =  "Final Flux-Gate Positions for the Fine-Tuned Descent"


tile_size = 0.262 #m
ntiles = [5, 9, 5]
totalsize = ntiles .* tile_size

g, vertex_positions = cuboid_system(totalsize, ntiles,
    skipfaces = [false, false, true, false, false, false])

g_2, vertex_positions_2 = cuboid_system([0.97, 0.97, 0.87], [6,6,6],
    skipfaces = [false, false, false, false, false, false])

# load data
poi  = readdlm("$(setup)/poi_$(file).csv" , ',' , Float64)[end-7:end,:]
poi_2 = readdlm("$(setup)/poi_$(file).csv" , ',' , Float64)[1:8,:]

# reformat final poi into similar format as poi in coils package
poi_vecs = []
poi_vecs_2 = []
for i in 0:7
    push!(poi_vecs, poi[end-i,:])
    push!(poi_vecs_2, poi_2[end-i,:])
end

# generate poi plot
plot_final(g, vertex_positions, poi_vecs, g_2, vertex_positions_2; standalone = true, alpha = 1, titles = titles)
gcf()
# save image
gcf().savefig("$(setup)/Images/$(titles)_initial_final.png", dpi = 1000)



## Initial and final positions of the algorithm
plot_initial_final(g, vertex_positions, poi_vecs, poi_vecs_2, g_2, vertex_positions_2; standalone = true, alpha = 1, titles = titles)
gcf()
# save image
gcf().savefig("$(setup)/Images/$(titles)_initial.png", dpi = 1000)


## How the positions change based on number of sensors

titles = "Final Flux-Gate Positions for varying Number of Gates"
final_positions = []
for i in 4:20
    final_positions_i = []
    for j in 0:i-1
        push!(final_positions_i,readdlm("$(setup)/num_sens_vals/poi_num_sens_$(i)_BB_5050.csv" , ',' , Float64)[end-j,:])
    end
    push!(final_positions, final_positions_i)
end
final_positions
plot_num_sens(g, vertex_positions, final_positions, g_2, vertex_positions_2; standalone = true, alpha = 1, titles = titles)
gcf()

# condition numbers after the num sens descent

## Tims Poi Descent
# Poi
file = "tims_tuned_SPSA_scaled"

titles =  "Final Flux-Gate Positions for the Fine-Tuned Descent"

tile_size = 0.262 #m
ntiles = [5, 9, 5]
totalsize = ntiles .* tile_size

g, vertex_positions = cuboid_system(totalsize, ntiles,
    skipfaces = [false, false, true, false, false, false])

g_2, vertex_positions_2 = cuboid_system([0.97, 0.97, 0.87], [6,6,6],
    skipfaces = [false, false, false, false, false, false])

# load data
poi  = readdlm("$(setup)/poi_$(file).csv" , ',' , Float64)[end-7:end,:]
poi_2 = readdlm("$(setup)/poi_$(file).csv" , ',' , Float64)[1:8,:]

# reformat final poi into similar format as poi in coils package
poi_vecs = []
poi_vecs_2 = []
for i in 0:7
    push!(poi_vecs, poi[end-i,:])
    push!(poi_vecs_2, poi_2[end-i,:])
end


plot_initial_final(g, vertex_positions, poi_vecs, poi_vecs_2, g_2, vertex_positions_2; standalone = true, alpha = 1, titles = titles)
gcf()
# save image
gcf().savefig("$(setup)/Images/$(titles)_initial.png", dpi = 1000)

# Conds

condition_num = readdlm("$(setup)/conds_$(file).csv" , ',' , Float64)

clf()

gca().plot(condition_num)
gcf().xtitle("Iteration (Arbitrary Units)")
gcf().ytitle("Condition Number (Arbitrary Units)")
gcf().title(titles)

gcf().savefig("$(setup)/Images/$(titles)_initial.png", dpi = 1000)

## Corner Poi Descent
# Poi

file = "tims_tuned_SPSA_scaled"

titles =  "Final Flux-Gate Positions for the Fine-Tuned Descent"

tile_size = 0.262 #m
ntiles = [5, 9, 5]
totalsize = ntiles .* tile_size

g, vertex_positions = cuboid_system(totalsize, ntiles,
    skipfaces = [false, false, true, false, false, false])

g_2, vertex_positions_2 = cuboid_system([0.97, 0.97, 0.87], [6,6,6],
    skipfaces = [false, false, false, false, false, false])

# load data
poi  = readdlm("$(setup)/poi_$(file).csv" , ',' , Float64)[end-7:end,:]
poi_2 = readdlm("$(setup)/poi_$(file).csv" , ',' , Float64)[1:8,:]

# reformat final poi into similar format as poi in coils package
poi_vecs = []
poi_vecs_2 = []
for i in 0:7
    push!(poi_vecs, poi[end-i,:])
    push!(poi_vecs_2, poi_2[end-i,:])
end


plot_initial_final(g, vertex_positions, poi_vecs, poi_vecs_2, g_2, vertex_positions_2; standalone = true, alpha = 1, titles = titles)
gcf()
# save image
gcf().savefig("$(setup)/Images/$(titles)_initial.png", dpi = 1000)

# Conds

condition_num = readdlm("$(setup)/conds_$(file).csv" , ',' , Float64)

clf()

gca().plot(condition_num)
gcf().xtitle("Iteration (Arbitrary Units)")
gcf().ytitle("Condition Number (Arbitrary Units)")
gcf().title(titles)

gcf().savefig("$(setup)/Images/$(titles)_initial.png", dpi = 1000)
## Side Poi Descent
# Poi

file = "tims_tuned_SPSA_scaled"

titles =  "Final Flux-Gate Positions for the Fine-Tuned Descent"

tile_size = 0.262 #m
ntiles = [5, 9, 5]
totalsize = ntiles .* tile_size

g, vertex_positions = cuboid_system(totalsize, ntiles,
    skipfaces = [false, false, true, false, false, false])

g_2, vertex_positions_2 = cuboid_system([0.97, 0.97, 0.87], [6,6,6],
    skipfaces = [false, false, false, false, false, false])

# load data
poi  = readdlm("$(setup)/poi_$(file).csv" , ',' , Float64)[end-7:end,:]
poi_2 = readdlm("$(setup)/poi_$(file).csv" , ',' , Float64)[1:8,:]

# reformat final poi into similar format as poi in coils package
poi_vecs = []
poi_vecs_2 = []
for i in 0:7
    push!(poi_vecs, poi[end-i,:])
    push!(poi_vecs_2, poi_2[end-i,:])
end


plot_initial_final(g, vertex_positions, poi_vecs, poi_vecs_2, g_2, vertex_positions_2; standalone = true, alpha = 1, titles = titles)
gcf()
# save image
gcf().savefig("$(setup)/Images/$(titles)_initial.png", dpi = 1000)

# Conds

condition_num = readdlm("$(setup)/conds_$(file).csv" , ',' , Float64)

clf()

gca().plot(condition_num)
gcf().xtitle("Iteration (Arbitrary Units)")
gcf().ytitle("Condition Number (Arbitrary Units)")
gcf().title(titles)

gcf().savefig("$(setup)/Images/$(titles)_initial.png", dpi = 1000)


## Fine tuned Poi Descent
# Poi

file = "fine_tuned_SPSA_scaled"

titles =  "Final Flux-Gate Positions for the Fine-Tuned Descent"

tile_size = 0.262 #m
ntiles = [5, 9, 5]
totalsize = ntiles .* tile_size

g, vertex_positions = cuboid_system(totalsize, ntiles,
    skipfaces = [false, false, true, false, false, false])

g_2, vertex_positions_2 = cuboid_system([0.97, 0.97, 0.87], [6,6,6],
    skipfaces = [false, false, false, false, false, false])

# load data
poi  = readdlm("$(setup)/poi_$(file).csv" , ',' , Float64)[end-7:end,:]
poi_2 = readdlm("$(setup)/poi_$(file).csv" , ',' , Float64)[1:8,:]

# reformat final poi into similar format as poi in coils package
poi_vecs = []
poi_vecs_2 = []
for i in 0:7
    push!(poi_vecs, poi[end-i,:])
    push!(poi_vecs_2, poi_2[end-i,:])
end


plot_initial_final(g, vertex_positions, poi_vecs, poi_vecs_2, g_2, vertex_positions_2; standalone = true, alpha = 1, titles = titles)
gcf()
# save image
gcf().savefig("$(setup)/Images/$(titles)_initial.png", dpi = 1000)

# Conds

condition_num = readdlm("$(setup)/conds_$(file).csv" , ',' , Float64)

clf()

gca().plot(condition_num)
gcf().xtitle("Iteration (Arbitrary Units)")
gcf().ytitle("Condition Number (Arbitrary Units)")
gcf().title(titles)

gcf().savefig("$(setup)/Images/$(titles)_initial.png", dpi = 1000)


## Condition number plots

## For the different number of sensors


clf()
condlist = []
for i in 4:20
    conds = readdlm("OutputBabySFC/COMSOL_Fields/num_sens_vals/conds_num_sens_$(i)_BB_5050.csv" , ',' , Float64)
    gca().plot(conds)
    push!(condlist, (last(conds), i))
end
condlist


gcf()

## For the AIRSS Endpoints

poi = readdlm("OutputBabySFC/COMSOL_Fields/airss_end_conds_5050.csv" , ',' , Float64)


conds = readdlm("OutputBabySFC/COMSOL_Fields/airss_end_conds_5050.csv" , ',' , Float64)

clf()
#plot(conds[:,1])
##
scatter((conds[:,1]), log.(conds[:,2]))

gcf()

clf()
bins = range(0,stop = 40, step = 1)
gcf().hist2D(log.(conds[:,1]), bins= bins, title = "Condition Numbers after AIRSS Descent") ## plot after 200 steps 100 initialisations~
histogram(log.(conds[:,2]), bins= bins, title = "Condition Numbers at the Start of AIRSS Initialisations") ## plot before 200 steps 100 initialisations

gcf()




## For the Fine Tuning


conds  = readdlm("OutputBabySFC/COMSOL_Fields/cond_fine_tuned_SPSA_5050.csv" , ',' , Float64)




## Field Stability
