current_path = pwd()
current_path = "/Users/LTV/dev_projects/innoPositioning-Julia"
push!(LOAD_PATH, "$(current_path)/src")

using MapPrimitives
using ImageTree
using MapPlan
using MapBuilder
using HDF5, JLD
using Plots

walls = MapPlan.read_data_2d("$(current_path)/res/coverage/walls2d.txt")

ssm = load("$(current_path)/res/coverage/full.jld","data")

plot(ssm,seriestype=:heatmap)
