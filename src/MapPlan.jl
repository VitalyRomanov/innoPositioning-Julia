module MapPlan
  using Geometry
  using MapPrimitives
  # using RadixTree
  using PBSM

  export mapPlan
  export read_data,plot_walls,downsample_maps,no_walls_on_path,walls_on_path

  type mapPlan
    walls::Array{Wall3D}
    # AP::Array{Float64}
    limits::Array{Int}
    # index::RadixTree.radixTree
    index::PBSM.Pbsm
    vis_matr::Array{Bool}
  end

  function query_walls(path::Line,wall_index)
    # return RadixTree.probe(wall_index,line2mbr(path))
    return PBSM.probe(wall_index,line2mbr(path))
  end

  function create_index(walls,lims,grid_scl = 10.)
    # return RadixTree.create_index(RadixTree.obj2mbr(walls,wall2mbr),lims)
    return PBSM.create_index(PBSM.obj2mbr(walls,wall2mbr),lims,grd_scl = grid_scl)
  end



  function downsample_maps(signal_maps,factor)
    println("Downsampling maps")
    new_maps = Array{Array{Float64}}(length(signal_maps))
    for (ind,map) in enumerate(signal_maps)
      new_maps[ind] = downsample_map(map,2)
    end
    return new_maps
  end

  function downsample_map(signal_map,factor)
    new_size = [Int(floor(size(signal_map,1)/factor)),Int(floor(size(signal_map,2)/factor))]
    new_map = Array{Float64}(new_size[1],new_size[2])
    for i=1:new_size[1]
      for j=1:new_size[2]
        ii = (i-1)*factor+1
        jj = (j-1)*factor+1
        new_map[i,j] = maximum(signal_map[ii:ii+factor-1,jj:jj+factor-1])
      end
    end
    return new_map
  end

  function plot_walls_2d(walls,x_max,y_max)
    gr()
    town_plan = plot([0,0],[0,0],linecolor = :black,xlims=(0,x_max),ylims=(0,y_max),grid=false,legend=false,axis=false)
    for wall in walls
      plot!(wall[:,1],wall[:,2],linecolor = :black,xlims=(0,x_max),ylims=(0,y_max),grid=false,legend=false,axis=false)
    end
    return town_plan
  end

  function plot_walls(walls,x_max,y_max)
    gr()
    town_plan = plot([0,0],[0,0],linecolor = :black,xlims=(0,x_max),ylims=(0,y_max),grid=false,legend=false,axis=false)
    for wall in walls
      plot!([wall.points[1].val[1],wall.points[4].val[1]],[wall.points[1].val[2],wall.points[4].val[2]],linecolor = :black,xlims=(0,x_max),ylims=(0,y_max),grid=false,legend=false,axis=false)
    end
    return town_plan
  end

  function read_data_2d(data_file)
    walls = []
    input = open(data_file)
    for (line_ind,line) in enumerate(eachline(input))
      strvec = split(line[1:end-2]," ")
      v = map(x->parse(Float64,x),strvec)
      wall = reshape(v[:],2,Int(length(v)/2))'
      push!(walls,wall)
    end
    return walls
  end


  function create_all_vertex_paths(polygon1::Array{Array{Float64}},polygon2::Array{Array{Float64}})
    paths = Array{Line}(length(polygon1)*length(polygon2))
    path_counter = 1
    for i = 1:length(polygon1)
      for j = 1:length(polygon2)
        paths[path_counter] = Line(polygon1[i],polygon2[j])
        path_counter += 1
      end
    end
    return paths
  end



  function shortest_polygon_distance(polygon1::Array{Array{Float64}},polygon2::Array{Array{Float64}})
    min_dist  = 1.e50

    for i = 1:length(polygon1)
      for j = 1:length(polygon2)
        min_dist = min(min_dist,norm(polygon1[i]-polygon2[j]))
      end
    end

    return min_dist
  end


  function create_ap_visibility(plan::mapPlan,AP)
    ap_visibility = Array{Bool}(length(plan.walls))*false
    ap_list = push!(Array{Array{Float64}}(0),AP)

    for (couple_ind,couple_wall) in enumerate(plan.walls)
      shrkd_wll_plgn2 = shrink_polygon(couple_wall.polygon,.2)

      if shortest_polygon_distance(ap_list,shrkd_wll_plgn2) > 200.
        continue
      end

      result = false
      test_path = Line()
      test_path.v1 = AP
      for vertex in shrkd_wll_plgn2
        test_path.v2 = vertex
        if no_walls_on_path(test_path,plan)
          result = true
          break
        end
      end
      # paths = create_all_vertex_paths(ap_list,shrkd_wll_plgn2)
      # result = false
      # for path in paths
      #   if no_walls_on_path(path,plan)
      #     result = true
      #     break
      #   end
      # end
      ap_visibility[couple_ind] = result
    end
    return ap_visibility
  end



  function create_wall_visibility_matrix(plan::mapPlan)
    function initfc(S::SharedArray)
        for s in eachindex(S)
            S[s] = false
        end
    end
    visibility_matrix = SharedArray{Bool}((length(plan.walls),length(plan.walls)),init=initfc)
    # visibility_matrix = Array(Bool,length(plan.walls),length(plan.walls))*false
    # for (wall_id,wall) in enumerate(plan.walls)
    @sync @parallel for wall_id = 1:length(plan.walls)
      print("\rInspecting wall $(wall_id)/$(length(plan.walls))...    ")
      wall = plan.walls[wall_id]
      shrkd_wll_plgn1 = shrink_polygon(wall.polygon,.2)
      for (couple_ind,couple_wall) in enumerate(plan.walls[wall_id+1:end])

        shrkd_wll_plgn2 = shrink_polygon(couple_wall.polygon,.2)

        if shortest_polygon_distance(shrkd_wll_plgn1,shrkd_wll_plgn2) > 200.
          continue
        end

        test_path = Line(); result = false
        for v1 in shrkd_wll_plgn1, v2 in shrkd_wll_plgn2
          test_path.v1 = v1; test_path.v2 = v2;
          if no_walls_on_path(test_path,plan)
            result = true
            break
          end
        end


        # paths = create_all_vertex_paths(shrkd_wll_plgn1,shrkd_wll_plgn2)
        #
        # result = false
        # for path in paths
        #   if no_walls_on_path(path,plan)
        #     result = true
        #     break
        #   end
        # end
        visibility_matrix[wall_id,wall_id+couple_ind] = result
      end
    end

    for i in 1:size(visibility_matrix,1)
      for j in i+1:size(visibility_matrix,1)
        visibility_matrix[j,i] = visibility_matrix[i,j]
      end
      visibility_matrix[i,i] = false
    end

    print("\n")

    return convert(Array,visibility_matrix)
  end


  function no_walls_on_path(path::Line,plan::mapPlan)
    shrink_line!(path,float_err_marg)
    filtered_walls = plan.walls[query_walls(path,plan.index)]
    for wall in filtered_walls
      if MapPrimitives.get_intersection_point(path,wall)!=-1
        return false
      end
    end
    return true
  end

  function walls_on_path(path::Line,plan::mapPlan)
    shrink_line!(path,float_err_marg)
    wop = 0
    # for wall_id in query_walls(path,plan.index)
    #     if MapPrimitives.get_intersection_point(path,plan.walls[wall_id])!=-1
    #         wop += 1
    #     end
    # end
    filtered_walls = plan.walls[query_walls(path,plan.index)]
    for wall in filtered_walls
      if MapPrimitives.get_intersection_point(path,wall)!=-1
        wop += 1
      end
    end
    return wop
  end

end
