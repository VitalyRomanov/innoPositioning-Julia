push!(LOAD_PATH, "/Users/LTV/Dropbox/mapGen/fast_vis_matr")
cd("/Users/LTV/Dropbox/mapGen/fast_vis_matr/2ref")
using HDF5,JLD
using MapPrimitives
using MapPlan
using MapBuilder
using PyPlot
using ImageFiltering

walls,bbox = MapPlan.read_data("../walls.txt")
signal_strength_matrix = load("full.jld","data")
for i in 1:size(signal_strength_matrix,1)
  for j in 1:size(signal_strength_matrix,2)
    if signal_strength_matrix[i,j] < -110.
      signal_strength_matrix[i,j] = -100.
    end
  end
end
blu = imfilter(signal_strength_matrix, Kernel.gaussian(5));


pcolormesh(blu')
clim(-110,-30)

for wall in walls
  plot([wall.points[1].val[1],wall.points[4].val[1]],[wall.points[1].val[2],wall.points[4].val[2]],"w")
end
