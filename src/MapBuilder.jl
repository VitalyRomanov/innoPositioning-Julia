module MapBuilder
  using MapPrimitives
  import MapPlan
  using ImageTree
  using Geometry

  export build_image_tree,calculate_signal_strength

  function create_wall_visibility_matrix(wall_ind::WallIndex)
    visibility_matrix = Array(Bool,length(wall_ind.walls),length(wall_ind.walls))
    tic()
    for (wall_index,wall) in enumerate(wall_ind.walls)
      toc()
      tic()
      print("Inspecting wall ",wall_index,"...\n")
      w1v1 = wall.points[1].val[1:2]
      w1v2 = wall.points[4].val[1:2]
      for (couple_ind,couple_wall) in enumerate(wall_ind.walls[wall_index+1:end])
        # print("\t checking friend ",couple_ind,"\n")
        w2v1 = couple_wall.points[1].val[1:2]
        w2v2 = couple_wall.points[4].val[1:2]

        paths = Array(Line,4)
        paths[1] = Line(w1v1,w2v1)
        paths[2] = Line(w1v1,w2v2)
        paths[3] = Line(w1v2,w2v1)
        paths[4] = Line(w1v2,w2v2)

        result = false
        for path in paths
          if no_walls_on_path(path,wall_ind)
            result = true
            break
          end
        end
        visibility_matrix[wall_index,wall_index-1+couple_ind] = result
      end
    end
    return visibility_matrix
  end


  function calculate_offsprings(image::treeNode,image_id::Int,tree_size::Int,wall_index::WallIndex,AP::Point,visibility_matrix;max_levels = 2, distance_threshold = 150)
    offsprings = Array(treeNode,length(wall_index.walls))
    children = Array(Int,length(wall_index.walls))
    real_offsprings = 0
    current_offset = 0

    if image.level < max_levels
      if image.assigned_wall == -1
        feasible_walls = wall_index.walls
      else
        feasible_walls = wall_index.walls[find(visibility_matrix[image.assigned_wall,:])]
      end
      for wall in feasible_walls
        if norm(wall.points[1].val-AP.val)<distance_threshold
          if wall.id != image.assigned_wall
            paths = Array(Line,2)
            paths[1] = Line(image.location.val[1:2],wall.points[1].val[1:2])
            paths[2] = Line(image.location.val[1:2],wall.points[4].val[1:2])

            wall_visible = false
            for path in paths
              if no_walls_on_path(path,wall_index)
                wall_visible = true
                break
              end
            end
            if wall_visible
              position = reflection_from_plane(image.location.val,wall.plane_equation)
              new_image = MapPrimitives.Point(position)
              new_node = ImageTree.treeNode(new_image,image.level+1,image_id,Array(Int,0),wall.id)
              real_offsprings += 1
              current_offset+=1
              offsprings[real_offsprings] = new_node
              children[real_offsprings] = tree_size+current_offset
            end
          end
        end
      end
    end
    # potentially there can be a problem of not all valuable nodes included in the tree
    image.children = children[1:real_offsprings]
    return offsprings[1:real_offsprings]
  end



  function build_image_tree(AP::Point,wall_index::WallIndex,visibility_matrix; pathloss_distance_threshold = 180)
    image_tree = Array(treeNode,0)
    root = treeNode(AP,0,-1,[],-1)
    # recalculate the image tree for every reception point to avoid combinatorial explosion
    push!(image_tree,root)
    tree_position = 1
    iteration_steps = length(image_tree)
    while tree_position <= iteration_steps
      offsprings = calculate_offsprings(image_tree[tree_position],tree_position,length(image_tree),wall_index,AP,visibility_matrix, distance_threshold = pathloss_distance_threshold)
      append!(image_tree,offsprings)
      tree_position += 1
      iteration_steps = length(image_tree)
    end
    print("Image tree is built: ",length(image_tree),"\n")
    return image_tree
  end



  function get_path_vertices(image_tree,node_ind,Rx::Point,wall_index::WallIndex)
    vertices = Array(Array{Float64},0)
    associated_walls = Array(Int,0)

    image = image_tree[node_ind]

    push!(vertices,Rx.val)
    push!(associated_walls,0)

    while image.parent!=-1
      current_image_loc = image.location.val
      current_wall = image.assigned_wall

      reflection_point = get_intersection_point(vertices[end],current_image_loc,wall_index.walls[current_wall])

      if reflection_point == -1
        return [],[],[]
      end

      if !no_walls_on_path(Line(vertices[end][1:2],reflection_point[1:2]),wall_index)
        return [],[],[]
      end

      push!(vertices,reflection_point)
      push!(associated_walls,current_wall)
      image = image_tree[image.parent]
    end
    if !no_walls_on_path(Line(vertices[end][1:2],image.location.val[1:2]),wall_index)
      return [],[],[]
    end
    push!(vertices,image.location.val)
    push!(associated_walls,image.assigned_wall)
    return vertices,associated_walls
  end


  function lines_crossed(line1::Line,line2::Line)
    a = line1.v1[2] - line1.v2[2]
    b = line1.v2[1] - line1.v1[1]
    c = (line1.v2[2] - line1.v1[2])*line1.v1[1] - (line1.v2[1] - line1.v1[1])*line1.v1[2]

    d = line2.v1[2] - line2.v2[2]
    e = line2.v2[1] - line2.v1[1]
    f = (line2.v2[2] - line2.v1[2])*line2.v1[1] - (line2.v2[1] - line2.v1[1])*line2.v1[2]

    den = (-d*b + a*e)

    if den == 0
      return false
    end

    y = (c*d - f*a) / den
    x = (b*f - c*e) / den

    ip = [x,y]

    if dot(line1.v1-ip,line1.v2-ip)<=0 && dot(line2.v1-ip,line2.v2-ip)<=0
      return true
    else
      return false
    end
  end

  function sector_crossed(ray::Line,sector::Sector)
    for edge in sector.geometry
      if lines_crossed(edge,ray)
        return true
      end
    end
    return false
  end



  function query_walls(v1,v2,wall_index)
    ray = Line(v1[1:2],v2[1:2])

    x1 = Int(ceil(v1[1]/wall_index.sector_size))
    x2 = Int(ceil(v2[1]/wall_index.sector_size))
    y1 = Int(ceil(v1[2]/wall_index.sector_size))
    y2 = Int(ceil(v2[2]/wall_index.sector_size))

    ind1 = (y1-1)*wall_index.sectors_dim[1]+x1
    ind2 = (y2-1)*wall_index.sectors_dim[1]+x2

    if x1<1 || y1<1 || x1 >= wall_index.sectors_dim[1] || y1 >= wall_index.sectors_dim[2]
     ind1 = -1
    end
    if x2<1 || y2<1 || x2 >= wall_index.sectors_dim[1] || y2 >= wall_index.sectors_dim[2]
     ind2 = -1
    end

    crossed_sectors = Set()

    if ind1 != -1
      union!(crossed_sectors,ind1)
    end
    if ind2 != -1
      union!(crossed_sectors,ind2)
    end

    if ind1!=ind2
      for (sector_ind,sector) in enumerate(wall_index.sectors)
        if sector_crossed(ray,sector)
          union!(crossed_sectors,sector_ind)
        end
      end
    end

    wi = Set()
    for sector_ind in collect(crossed_sectors)
      union!(wi,wall_index.sectors[sector_ind].walls)
    end

    return collect(wi)
  end


  function no_walls_on_path(path::Line,wall_index::WallIndex)
    filtered_walls = wall_index.walls[query_walls(path.v1,path.v2,wall_index)]
    for wall in filtered_walls
      if lines_crossed(path,Line(wall.points[1].val[1:2],wall.points[4].val[1:2]))
        return false
      end
    end
    return true
  end


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


  function calculate_signal_strength(Rx::Point,image_tree::Array{treeNode},wall_index::WallIndex; attenuation_exponent = 2.2,reflection_coef = 12.53,transmission_coef = 100., P0 = -30.)
    distance = Array(Float64,0)
    path_loss = Array(Float64,0)
    for (node_ind,image) in enumerate(image_tree)
      path_vertices,associated_walls = get_path_vertices(image_tree,node_ind,Rx,wall_index)

      if length(path_vertices)>0
        dist_sum = get_path_distance(path_vertices)

        pl = 10^((P0-10*attenuation_exponent*log10(dist_sum)-reflection_coef*(length(associated_walls)-2))/10)

        push!(distance,dist_sum)
        push!(path_loss,pl)
      end
    end

    if sum(path_loss)==0
      return -900
    else
      return Int(round(10*log10(sum(path_loss))))
    end
  end

  function get_path_distance(path)
    dist = 0.
    for i in 1:(length(path)-1)
      dist += norm(path[i]-path[i+1])
    end
    return dist
  end


end
