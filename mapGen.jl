# plane equation does not look the same, but the reflections are calculated correctly -- resolved
# why do you need check_if_point_belongs_to_line_segment(p1, p2, p)? -- resolved, to chech whether it belongs to the path
# check_if_point_belongs_to_wall checks only two walls? can be improved
# signal strength increases when transmissions are introduced
# optimize math
# in order to make everything parallel I need to transfer data to all the workers

current_path = pwd()
push!(LOAD_PATH, "$(current_path)/src")
# push!(LOAD_PATH, "/Users/LTV/Dropbox/mapGen/fast_vis_matr")

using MapPrimitives
using ImageTree
using MapPlan
using MapBuilder
using HDF5, JLD
using PyPlot






function calculate_signal_strength_matrix(image_tree,wall_index,x_max,y_max,AP_ind)
  print("Calculating signal strength matrix\n")
  # signal_strength_matrix = ones(Float64,x_max,y_max)*-1000
  ss = ones(Float64,x_max,y_max)*-900
  # paths = Array(Array{Float64},x_max,y_max)

  pathloss_distance_threshold = 150.

  path_dump = open("$(current_path)/res/coverage/path_dump_$(AP_ind).txt","w")

  for x in collect(1:x_max)
    # tic()
    for y in collect(1:y_max)
      if norm([x,y,1]-image_tree[1].location.val) < pathloss_distance_threshold
        # signal_strength_matrix[x,y] =  MapBuilder.calculate_signal_strength(MapPrimitives.Point([x,y,1.0]),image_tree,wall_index)
        ss[x,y],pp =  MapBuilder.calculate_signal_strength(MapPrimitives.Point([x,y,1.0]),image_tree,wall_index)
        if length(pp)>0
          write(path_dump,"$(x) $(y) $(pp)\n")
        end
      end
    end
    # toc()
    print("\r$(image_tree[1].location.val[1:2]) $(x)/$(x_max)                    ")
  end
  close(path_dump)
  print("\n")
  return ss
end

# AP = MapPrimitives.Point([599.344866581,466.28583739,3.0])



APs = [MapPrimitives.Point([200.99503661, 394.493304064,3.0]),
      MapPrimitives.Point([217.939597303, 372.489304064,3.0]),
      MapPrimitives.Point([377.462713731, 443.725770731,3.0]),
      MapPrimitives.Point([366.624301036, 429.500704064,3.0]),
      MapPrimitives.Point([597.344866581, 468.285837397,3.0]),
      MapPrimitives.Point([523.170670315, 337.260037397,3.0]),
      MapPrimitives.Point([836.721369973, 115.771304064,3.0]),
      MapPrimitives.Point([374.409639732, 243.908037397,3.0]),
      MapPrimitives.Point([615.831466173, 384.380570731,3.0]),
      MapPrimitives.Point([625.372322419, 511.628237397,3.0])]

# AP = MapPrimitives.Point([377.462713731,443.725770731,3.0])

ssms = []

walls,bbox = MapPlan.read_data("$(current_path)/res/coverage/walls.txt")
x_max = Int64(round(bbox[2].points[3].val[1]))
y_max = Int64(round(bbox[2].points[3].val[2]))
wall_index = MapPrimitives.create_wall_index(walls,x_max,y_max)
# visibility_matrix = MapBuilder.create_wall_visibility_matrix(wall_index)
# save("vis_ind.jld", "visibility_matrix", visibility_matrix)
visibility_matrixx = load("$(current_path)/res/coverage/vis_ind.jld")
visibility_matrix = visibility_matrixx["visibility_matrix"]

for i in 1:size(visibility_matrix,1)
  for j in i+1:size(visibility_matrix,1)
    visibility_matrix[j,i] = visibility_matrix[i,j]
  end
  visibility_matrix[i,i] = false
end

for (AP_ind,AP) in enumerate(APs)
  x_max = Int64(round(bbox[2].points[3].val[1]))
  y_max = Int64(round(bbox[2].points[3].val[2]))

  wall_index = MapPrimitives.create_wall_index(walls,x_max,y_max)
  image_tree = MapBuilder.build_image_tree(AP,wall_index,visibility_matrix,pathloss_distance_threshold = 150)
  signal_strength_matrix = calculate_signal_strength_matrix(image_tree,wall_index,x_max,y_max,AP_ind)
  pcolormesh(signal_strength_matrix)
  clim(-120,-30)
  # colorbar()
  savefig("$(AP_ind).png")
  save("$(current_path)/res/coverage/data_$(AP_ind).jld", "data", signal_strength_matrix)

  push!(ssms,signal_strength_matrix)
end

signal_strength_matrix = ones(Float64,size(ssms[1],1),size(ssms[1],2))*-1000
for x in 1:size(ssms[1],1)
  for y in 1:size(ssms[1],2)
    for i in 1:length(ssms)
      signal_strength_matrix[x,y] = max(ssms[i][x,y],signal_strength_matrix[x,y])
    end
  end
end
pcolormesh(signal_strength_matrix)
clim(-120,-30)
# colorbar()
savefig("$(current_path)/res/coverage/full.png")
save("$(current_path)/res/coverage/full.jld", "data", signal_strength_matrix)
