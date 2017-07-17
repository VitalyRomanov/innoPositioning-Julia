current_path = pwd()
current_path = "/Users/LTV/dev_projects/innoPositioning-Julia"
push!(LOAD_PATH, "$(current_path)/src")

using MapPrimitives
# using ImageTree
using MapPlan
# using MapBuilder
using HDF5, JLD
using Plots
# using PyPlot

walls,x_lims,y_lims = MapPlan.read_data("$(current_path)/res/coverage/walls_tables.txt")

ssm = load("$(current_path)/res/coverage/full.jld","data")'

gr()

op = plot(ssm,seriestype=:heatmap,seriescolor=ColorGradient([colorant"white", colorant"orange", colorant"red"]),zlims=(-40,30),legend = false,grid=false,axis=false)

for wall in walls
  temp_x = [wall.points[1].val[1],wall.points[4].val[1]]
  temp_y = [wall.points[1].val[2],wall.points[4].val[2]]*2
  plot!(temp_x,temp_y,linecolor=:black,xlims=(x_lims[1],x_lims[2]),ylims=(x_lims[1],x_lims[2]))
end

# for wall in walls
#   plot!(wall[:,1],wall[:,2],linecolor=:black,xlims=(0,1038),ylims=(0,627))
# end

plot!(op)
# savefig("coverage_faded.png")
