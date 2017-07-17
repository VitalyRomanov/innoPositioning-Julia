# current_path = pwd()
current_path = "/Users/LTV/dev_projects/innoPositioning-Julia"
push!(LOAD_PATH, "$(current_path)/src")

using MapPrimitives
using ImageTree
using MapPlan
using MapBuilder
using HDF5, JLD
using Plots
using RadixTree


# rx = [35.,17.]
tx = [30.,10.]
rxs = [
      [32.6,10.5],
      [34.8,9.4],
      [32.6,8.2],
      [32.6,5.9],
      [32.6,3.6],
      [29.2,5.9],
      [29.2,3.6],
      [20.,5.],
      [20.,1.],
      [14.,3.],
      [8.,3.],
      [35.,16.]]


function wall2mbr(wall::Wall)
  v1 = minimum([wall.points[1].val wall.points[2].val wall.points[3].val wall.points[4].val],2)
  v2 = maximum([wall.points[1].val wall.points[2].val wall.points[3].val wall.points[4].val],2)
  return MBR(v1,v2)
end


function calculate_signal_strength_matrix(image_tree,plan,rx;pl_dist_thr = 100.)
  # print("Calculating signal strength matrix\n")
  ss = Array(Float64,plan.limits[1,2],plan.limits[2,2])
  pv = []
  for x in [Int(round(rx[1]))]#collect(1:x_max)
    for y in [Int(round(rx[2]))]#collect(1:y_max)
      if norm([x,y,1]-image_tree[1].location.val) < pl_dist_thr
        ss[x,y],~,pv =  MapBuilder.calculate_signal_strength(MapPrimitives.Point([x,y,1.0]),image_tree,plan)
      end
    end
  end
  print("\n")
  return ss,pv
end



AP = MapPrimitives.Point([tx[1], tx[2],1.0])

ssms = []

walls,x_lims,y_lims,lims = MapPlan.read_data("$(current_path)/res/coverage/walls_tables.txt")

# lims = [x_lims';y_lims';[0,6]']

# wall_index = MapPrimitives.create_wall_index(walls,x_lims[2],y_lims[2])
wirt = RadixTree.create_radix_tree(RadixTree.obj2mbr(walls,wall2mbr),lims)

plan = MapPlan.mapPlan(walls,AP,lims,wirt,Array(Bool))


visibility_matrix = get_visibility_matrix("$(current_path)/res/coverage/vis_ind.jld",plan)

plan.vis_matr = visibility_matrix


image_tree = build_image_tree(plan)


for (rx_ind,rx) in enumerate(rxs)
  ss,pv = calculate_signal_strength_matrix(image_tree,plan,rx)

  paths = Array(Array{Float64},length(pv))

  for (p_ind,p) in enumerate(pv)
    paths[p_ind] = Array(Float64,length(p),2)
    for (v_ind,v) in enumerate(p)
      paths[p_ind][v_ind,1] = v[1]
      paths[p_ind][v_ind,2] = v[2]
    end
  end

  gr()

  town_plan = plot([0,0],[0,0],linecolor = :black,xlims=(x_lims[1],x_lims[2]),ylims=(y_lims[1],y_lims[2]),grid = false,legend = false, axis = false)

  for wall in walls
    temp_x = [wall.points[1].val[1],wall.points[4].val[1]]
    temp_y = [wall.points[1].val[2],wall.points[4].val[2]]
    # println("$(temp_x) $(temp_y)")
    plot!(temp_x,temp_y,linecolor=:black,xlims=(x_lims[1],x_lims[2]),ylims=(x_lims[1],x_lims[2]))
  end

  plot!([rx[1]],[rx[2]],markershape=:circle,markercolor=:red,annotations=(rx[1]-10,rx[2]+2,text("Rx",:left,:red,9,"Courier Bold")),xlims=(x_lims[1],x_lims[2]),ylims=(y_lims[1],y_lims[2]))
  plot!([tx[1]],[tx[2]],markershape=:circle,markercolor=:red,annotations=(tx[1]+5,tx[2]+2,text("Tx",:left,:red,9,"Courier Bold")),xlims=(x_lims[1],x_lims[2]),ylims=(y_lims[1],y_lims[2]))

  for path in paths
    plot!(path[:,1],path[:,2],linecolor = :red,xlims=(x_lims[1],x_lims[2]),ylims=(y_lims[1],y_lims[2]))
  end

  plot(town_plan,xlims=(x_lims[1],x_lims[2]),ylims=(y_lims[1],y_lims[2]),annotations=(10,10,text("$(ss[Int(round(rx[1])),Int(round(rx[2]))])",:left)))
  # println(ss[Int(round(rx[1])),Int(round(rx[2]))])
  println(norm(rx-tx))
  savefig("paths$(rx_ind).png")
end
