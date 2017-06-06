push!(LOAD_PATH, "/Users/LTV/Dropbox/mapGen/fast_vis_matr")
cd("/Users/LTV/Dropbox/mapGen/fast_vis_matr/")
using HDF5,JLD
using MapPrimitives
using MapPlan
using MapBuilder
using PyPlot



visibility_matrix = load("vis_ind.jld")
visibility_matrixx = visibility_matrix["s"]

for i in 1:size(visibility_matrixx,1)
  for j in i+1:size(visibility_matrixx,1)
    visibility_matrixx[j,i] = visibility_matrixx[i,j]
  end
end

walls,bbox = MapPlan.read_data("walls.txt")
x_max = Int64(round(bbox[2].points[3].val[1]))
y_max = Int64(round(bbox[2].points[3].val[2]))


for i in 1:100
  for wall in walls[find(visibility_matrixx[i,:])]
    plot([wall.points[1].val[1],wall.points[4].val[1]],[wall.points[1].val[2],wall.points[4].val[2]],"b")
  end
  plot([walls[i].points[1].val[1],walls[i].points[4].val[1]],[walls[i].points[1].val[2],walls[i].points[4].val[2]],"r")
  xlim([0,1000])
  ylim([0,600])
  savefig("wall_accosiations/$(i).png")
  close()
end

wall_index = MapPrimitives.create_wall_index(walls,x_max,y_max,100.)
#
# for wall in walls[wall_index.sectors[3].walls]
#   plot([wall.points[1].val[1],wall.points[4].val[1]],[wall.points[1].val[2],wall.points[4].val[2]])
# end

ray1 = MapPrimitives.Line([15,15],[50,15])
# ray2 = MapPrimitives.Line([15,15],[150,15])
ray3 = MapPrimitives.Line([15,15],[150,315])
#
qw = @time MapBuilder.query_walls(ray3.v1,ray3.v2,wall_index)
