# plane equation does not look the same, but the reflections are calculated correctly -- resolved
# why do you need check_if_point_belongs_to_line_segment(p1, p2, p)? -- resolved, to chech whether it belongs to the path
# check_if_point_belongs_to_wall checks only two walls? can be improved
# signal strength increases when transmissions are introduced
# optimize math
# in order to make everything parallel I need to transfer data to all the workers

current_path = pwd()
push!(LOAD_PATH, "$(current_path)/src")

using MapPrimitives
using ImageTree
using RadixTree
using MapPlan
using MapBuilder
using HDF5, JLD
using PyPlot




function calculate_signal_strength_matrix(image_tree,plan,AP_ind)
  println("Calculating signal strength matrix")
  # signal_strength_matrix = ones(Float64,x_max,y_max)*-1000
  ss = ones(Float64,plan.limits[1,2],plan.limits[2,2])*-900
  # paths = Array(Array{Float64},x_max,y_max)

  pathloss_distance_threshold = 150.

  path_dump = open("$(current_path)/res/coverage/path_dump_$(AP_ind).txt","w")

  for x in collect(1:plan.limits[1,2])
    for y in collect(1:plan.limits[2,2])
      if norm([x,y,1]-image_tree[1].location.val) < pathloss_distance_threshold
        ss[x,y],pp,~ =  MapBuilder.calculate_signal_strength(MapPrimitives.Point([x,y,1.0]),image_tree,plan)
        if length(pp)>0
          write(path_dump,"$(x) $(y) $(pp)\n")
        end
      end
    end
    # toc()
    print("\r$(image_tree[1].location.val[1:2]) $(x)/$(plan.limits[1,2])                    ")
  end
  close(path_dump)
  print("\n")
  return ss
end

# AP = MapPrimitives.Point([599.344866581,466.28583739,3.0])



# APs = [MapPrimitives.Point([200.99503661, 394.493304064,3.0]),
#       MapPrimitives.Point([217.939597303, 372.489304064,3.0]),
#       MapPrimitives.Point([377.462713731, 443.725770731,3.0]),
#       MapPrimitives.Point([366.624301036, 429.500704064,3.0]),
#       MapPrimitives.Point([597.344866581, 468.285837397,3.0]),
#       MapPrimitives.Point([523.170670315, 337.260037397,3.0]),
#       MapPrimitives.Point([836.721369973, 115.771304064,3.0]),
#       MapPrimitives.Point([374.409639732, 243.908037397,3.0]),
#       MapPrimitives.Point([615.831466173, 384.380570731,3.0]),
#       MapPrimitives.Point([625.372322419, 511.628237397,3.0])]

APs = [MapPrimitives.Point([30.,10.,1.0])]

ssms = []

walls,~,~,lims = MapPlan.read_data("$(current_path)/res/coverage/walls_tables.txt")

# wall_index = MapPrimitives.create_wall_index(walls,x_lims[2],y_lims[2])
wirt = RadixTree.create_radix_tree(RadixTree.obj2mbr(walls,wall2mbr),lims)

plan = MapPlan.mapPlan(walls,APs[1],lims,wirt,Array(Bool))

visibility_matrix = get_visibility_matrix("$(current_path)/res/coverage/vis_ind.jld",plan)
# visibility_matrix_path = "$(current_path)/res/coverage/vis_ind.jld"
# if isfile(visibility_matrix_path)
#   println("Loading visibility matrix")
#   visibility_matrix = load(visibility_matrix_path,"visibility_matrix")
# else
#   println("Visibility matrix not found. Calculating...")
#   visibility_matrix = MapBuilder.create_wall_visibility_matrix(wall_index)
#   save(visibility_matrix_path,"visibility_matrix",visibility_matrix)
# end


plan.vis_matr = visibility_matrix


for (AP_ind,AP) in enumerate(APs)

  # wall_index = MapPrimitives.create_wall_index(walls,x_lims[2],y_lims[2])
  # image_tree = MapBuilder.build_image_tree(AP,wall_index,visibility_matrix,pathloss_distance_threshold = 150)

  plan.AP = AP

  image_tree = build_image_tree(plan)
  signal_strength_matrix = calculate_signal_strength_matrix(image_tree,plan,AP_ind)
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
