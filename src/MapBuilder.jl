module MapBuilder
  using MapPrimitives
  using MapPlan
  using ImageTree
  using Geometry
  using RadixTree


  export build_image_tree,calculate_signal_strength,get_visibility_matrix




  # function get_visibility_matrix(visibility_matrix_path,plan::mapPlan)
  #   if isfile(visibility_matrix_path)
  #     println("Loading visibility matrix")
  #     visibility_matrix = load(visibility_matrix_path,"visibility_matrix")
  #   else
  #     println("Visibility matrix not found. Calculating...")
  #     visibility_matrix = create_wall_visibility_matrix(plan)
  #     save(visibility_matrix_path,"visibility_matrix",visibility_matrix)
  #   end
  #   return visibility_matrix
  # end







  # function wall_position_2D(test_pos::Point,reference_wall::Wall)
  #   slope = (reference_wall.points[4].val[2]-reference_wall.points[1].val[2])/(reference_wall.points[4].val[1]-reference_wall.points[1].val[1])
  #   bias = reference_wall.points[4].val[2] - slope * reference_wall.points[4].val[1]
  #   return test_pos.val[2] > slope * test_pos.val[2] + bias
  # end


  # function calculate_offsprings(image_tree::Array{treeNode},image_id::Int,plan::mapPlan;max_levels = 2, distance_threshold = 150)
  #   # This function calculates the offsprings of a particular image in the image tree.
  #   # It is crucial to reduce the number of irrelevant images as much a possible. For this reason several optimizations are done.
  #   # 1. If the node is the root node, the walls are checked for being visible from the standpoint of the access point.
  #   # 2. Otherwise visibility matrix is used to determine whether the wall is visible from the standpoint of reflection wall
  #   # 3. The images should be built only for walls that are on the same side as the previous secondary source
  #   tree_size = length(image_tree)
  #   offsprings = Array(treeNode,length(plan.walls))
  #   children = Array(Int,length(plan.walls))
  #   real_offsprings = 0
  #   current_offset = 0
  #
  #   image = image_tree[image_id]
  #
  #   if image.level < max_levels
  #     # if the image is the root image no prior knowledge of wall visibility is available
  #     if image.assigned_wall == -1
  #       feasible_walls = plan.walls
  #     else
  #       feasible_walls = plan.walls[find(plan.vis_matr[image.assigned_wall,:])]
  #
  #       parent_wall = image_tree[image.parent].assigned_wall
  #       if parent_wall == -1
  #         parent_reference = plan.AP
  #       else
  #         parent_reference = plan.walls[parent_wall].points[1]
  #       end
  #       current_wall = plan.walls[image.assigned_wall]
  #       parent_position = wall_position_2D(parent_reference,current_wall)
  #     end
  #
  #     # we are going to reduce the number of feasible walls
  #     temp_feasible_walls = Array(Wall,length(feasible_walls))
  #     temp_position = 1
  #
  #     for wall in feasible_walls
  #       if norm(wall.points[1].val-plan.AP.val)<distance_threshold
  #         if image.assigned_wall == -1
  #           truly_feasible =  no_walls_on_path(Line(image.location.val,(wall.points[1].val+wall.points[4].val)/2),plan)
  #         else
  #           truly_feasible = (parent_position==wall_position_2D(wall.points[1],current_wall))
  #         end
  #
  #         if truly_feasible
  #           temp_feasible_walls[temp_position] = wall
  #           temp_position += 1
  #         end
  #       end
  #     end
  #
  #     feasible_walls = temp_feasible_walls[1:temp_position-1]
  #
  #     for wall in feasible_walls
  #       position = reflection_from_plane(image.location.val,wall.plane_equation)
  #       new_image = MapPrimitives.Point(position)
  #       new_node = ImageTree.treeNode(new_image,image.level+1,image_id,Array(Int,0),wall.id)
  #       real_offsprings += 1
  #       current_offset += 1
  #       offsprings[real_offsprings] = new_node
  #       children[real_offsprings] = tree_size+current_offset
  #     end
  #   end
  #
  #   image.children = children[1:real_offsprings]
  #   return offsprings[1:real_offsprings]
  # end
  #
  #
  #
  # function build_image_tree(plan::mapPlan; pl_thr_dist = 180)
  #   image_tree = Array(treeNode,0)
  #
  #   root = treeNode(plan.AP,0,-1,[],-1)
  #   push!(image_tree,root)
  #
  #   tree_position = 1
  #   iteration_steps = length(image_tree)
  #
  #   while tree_position <= iteration_steps
  #
  #     offsprings = calculate_offsprings(image_tree,tree_position,plan, distance_threshold = pl_thr_dist)
  #     append!(image_tree,offsprings)
  #
  #     tree_position += 1
  #     iteration_steps = length(image_tree)
  #   end
  #
  #   print("Image tree is built: ",length(image_tree),"\n")
  #   return image_tree
  # end

  function get_path_info(image_tree::Array{treeNode},node_ind::Int,Rx::Array{Float64},plan::mapPlan)

    image = image_tree[node_ind]

    last_vertex = Rx
    wop = 0::Int
    nrw = image.level
    dist = 0.::Float64

    while image.parent!=-1
      current_image_loc = image.location
      current_wall = image.assigned_wall

      susp_path = Line(last_vertex,current_image_loc)

      reflection_point = get_intersection_point(susp_path,plan.walls[current_wall])

      if reflection_point != -1
        susp_path = Line(last_vertex,reflection_point)
        wop += walls_on_path(susp_path,plan)
        dist += norm(last_vertex - reflection_point)
      else
        return -1.,nrw,wop
      end

      last_vertex = reflection_point
      image = image_tree[image.parent]
    end

    susp_path = Line(last_vertex,image.location)
    wop += walls_on_path(susp_path,plan)
    dist += norm(last_vertex - image.location)

    return dist,nrw,wop
  end



  function get_path_vertices(image_tree,node_ind,Rx::Array{Float64},plan::mapPlan)

    image = image_tree[node_ind]

    vertices = Array(Array{Float64},image.level+2)
    associated_walls = Array(Int,image.level)
    # vertices = Array(Array{Float64},0)
    # associated_walls = Array(Int,0)

    vertices[end] = Rx
    # push!(vertices,Rx.val)
    # push!(associated_walls,0)
    vert_pos = length(vertices)-1
    assw_pos = length(associated_walls)

    while image.parent!=-1
      current_image_loc = image.location
      current_wall = image.assigned_wall

      susp_path = Line(vertices[vert_pos+1],current_image_loc)

      reflection_point = get_intersection_point(susp_path,plan.walls[current_wall])

      if reflection_point != -1
        susp_path = Line(vertices[vert_pos+1],reflection_point,plan)
        wop = walls_on_path(susp_path,plan)
        return [],[],[]
      end

      # if !no_walls_on_path(Line(vertices[end],reflection_point),plan)
      #   return [],[],[]
      # end

      push!(vertices,reflection_point)
      push!(associated_walls,current_wall)
      image = image_tree[image.parent]
    end

    if !no_walls_on_path(Line(vertices[end],image.location.val),plan)
      return [],[],[]
    end
    push!(vertices,image.location.val)
    push!(associated_walls,image.assigned_wall)
    return vertices,associated_walls
  end




  # function sector_crossed(ray::Line,sector::Sector)
  #   for edge in sector.geometry
  #     if lines_crossed(edge,ray)
  #       return true
  #     end
  #   end
  #   return false
  # end


  # function query_walls(path::Line,wall_index::RadixTree.radixTree)
  #   return RadixTree.probe(wall_index,line2mbr(path))
  # end


  # function query_walls(ray::Line,wall_index::WallIndex)
  #
  #   # vv1 = convert(Array{Int},ceil(minimum([ray.v1 ray.v2],2)))
  #   # vv2 = convert(Array{Int},ceil(maximum([ray.v1 ray.v2],2)))
  #
  #   ##### try
  #   vv = convert(Array{Int},ceil([minimum([ray.v1 ray.v2],2) maximum([ray.v1 ray.v2],2)]))
  #   ind = (vv - [0 0;1 1])' * [1., wall_index.sectors_dim[1]]
  #   #####
  #
  #   x1 = Int(ceil(vv1[1]/wall_index.sector_size))
  #   x2 = Int(ceil(vv2[1]/wall_index.sector_size))
  #   y1 = Int(ceil(vv1[2]/wall_index.sector_size))
  #   y2 = Int(ceil(vv2[2]/wall_index.sector_size))
  #
  #   ind1 = (y1-1)*wall_index.sectors_dim[1]+x1
  #   ind2 = (y2-1)*wall_index.sectors_dim[1]+x2
  #
  #   if x1<1 || y1<1 || x1 >= wall_index.sectors_dim[1] || y1 >= wall_index.sectors_dim[2]
  #    ind1 = -1
  #   end
  #   if x2<1 || y2<1 || x2 >= wall_index.sectors_dim[1] || y2 >= wall_index.sectors_dim[2]
  #    ind2 = -1
  #   end
  #
  #   crossed_sectors = Set()
  #
  #   if ind1 != -1
  #     union!(crossed_sectors,ind1)
  #   end
  #   if ind2 != -1
  #     union!(crossed_sectors,ind2)
  #   end
  #
  #   if ind1!=ind2
  #     for (sector_ind,sector) in enumerate(wall_index.sectors)
  #       if sector_crossed(ray,sector)
  #         union!(crossed_sectors,sector_ind)
  #       end
  #     end
  #   end
  #
  #   wi = Set()
  #   for sector_ind in collect(crossed_sectors)
  #     union!(wi,wall_index.sectors[sector_ind].walls)
  #   end
  #
  #   return collect(wi)
  # end






  # function get_intersected_walls(v1,v2,wall_index::WallIndex)
  #   iw = Array(Int,0)
  #   filtered_walls = wall_index.walls[query_walls(v1,v2,wall_index)]
  #   for (wall_ind,wall) in enumerate(filtered_walls)
  #     ip = get_intersection_point(v1,v2,wall)
  #     if ip != -1
  #       push!(iw,wall_ind)
  #       # return 1
  #     end
  #   end
  #   # print(length(iw)," walls intersected\n")
  #   return iw
  # end


  function caclulate_signal_strength_matrix(ssm::Array{Float64},image_tree::Array{treeNode},plan::mapPlan)
    # add support for multiple resolutions
    for x = 1:size(ssm,1)
      for y = 1:size(ssm,2)
        # incorporate cell size of type Float64
        x_loc = Float64(x + plan.limits[1,1]) # multiply by cell size
        y_loc = Float64(y + plan.limits[2,1]) # multiply by cell size
        ssm[x,y] = calculate_signal_strength([x_loc,y_loc,1],image_tree,plan)
      end
    end
    return ssm
  end


  function calculate_signal_strength(Rx::Array{Float64},image_tree::Array{treeNode},plan::mapPlan; ae = 2,rc = 12.53,tc = 5., P0 = 0., sc = 0.5, pldthr = 200.)

    distance = Array(Float64,0)
    reflections = Array(Int,0)
    intersections = Array(Float64,0)
    path_loss = Array(Float64,0)

    used_images = 0
    # all_paths = Array(Float64,0)

    # pruned_branches = Set()

    # pv = []
    for (node_ind,image) in enumerate(image_tree)
      if norm(image.location - Rx) > pldthr
        continue
      end

      dist,nrw,wop = get_path_info(image_tree,node_ind,Rx,plan)

      if dist!=-1.
        used_images += 1
        push!(distance,dist)
        push!(reflections,nrw)
        push!(intersections,wop)

        pl = 10^( ( P0 - (10*ae*log10(dist)+20*log10(2.4e9) - 147.55 + sc*rc*nrw + tc*wop)) /10.)
        push!(path_loss,pl)
      end

      # if length(path_vertices)>0
      #   dist_sum = get_path_distance(path_vertices)
      #   append!(all_paths,[dist_sum;length(path_vertices)])
      #   push!(pv,path_vertices)
      #
      #   # pl = 10^((P0-10*attenuation_exponent*log10(dist_sum/.1)-reflection_coef*(length(associated_walls)-2)+147.55-20*log10(2.4e9))/10)
      #   pl = 10^( ( P0 - (10*ae*log10(dist_sum)+20*log10(2.4e9) - 147.55 + sc*rc*(length(associated_walls)-2))) /10.)
      #   # pl = 10^( ( P0 - (20*log10(dist_sum)+20*log10(2.4e9) - 147.55) ) /10.)
      #
      #   push!(distance,dist_sum)
      #   push!(path_loss,pl)
      # end
    end

    # println("$(used_images/length(image_tree)) images used")

    pl_sum = sum(path_loss)
    if pl_sum==0
      value = -900.
    else
      value = 10*log10(pl_sum)
    end

    return value
  end

  function get_path_distance(path)
    dist = 0.
    for i in 1:(length(path)-1)
      dist += norm(path[i]-path[i+1])
    end
    return dist
  end


end
