module MapPlan
  using Geometry
  using MapPrimitives
  using PBSM

  export mapPlan
  export read_data,downsample_maps,no_walls_on_path,walls_on_path

  mutable struct mapPlan
    walls::Array{Wall3D}
    limits::Array{Int}
    index::PBSM.Pbsm
    vis_matr::Array{Bool}
  end


  function query_walls(path::Line,wall_index)
    return PBSM.probe(wall_index,line2mbr(path))
  end

  function create_index(walls,lims;grid_scl = 10.)
    return PBSM.create_index(PBSM.obj2mbr(walls,wall2mbr),lims,grd_scl = grid_scl)
  end


  # function visualizeWallVis(project, range = [])
  #     if range == []
  #         range = 1:length(plan.walls)
  #     end
  #
  #     wallVisOutPath = project.path_save_data+"/wall_Vis/"
  #
  #     if !isdir(wallVisOutPath)
  #         mkdir(wallVisOutPath)
  #     end
  #
  #     plot()
  #
  #     for wall_id = range
  #         plot_walls!(project)
  #         plot_walls!(project,
  #                   ids = find(project.plan.vis_matr[wall_id,:]),
  #                   col = :red)
  #         plot_walls!(project,
  #                   ids = [wall_id]),
  #                   col = :green)
  #         savefig(wallVisOutPath+"wall_$(wall_id).svg")
  #     end
  # end



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

  # function plot_walls_2d(walls,x_max,y_max)
  #   gr()
  #   town_plan = plot([0,0],[0,0],linecolor = :black,xlims=(0,x_max),ylims=(0,y_max),grid=false,legend=false,axis=false)
  #   for wall in walls
  #     plot!(wall[:,1],wall[:,2],linecolor = :black,xlims=(0,x_max),ylims=(0,y_max),grid=false,legend=false,axis=false)
  #   end
  #   return town_plan
  # end

  # function plot_walls(walls,x_max,y_max;ids = [],col = :black)
  #   gr()
  #   if ids = []
  #       ids = 1:length(walls)
  #   end
  #   town_plan = plot([0,0],[0,0],linecolor = col,xlims=(0,x_max),ylims=(0,y_max),grid=false,legend=false,axis=false)
  #   for wall in walls[ids]
  #     plot!([wall.points[1].val[1],wall.points[4].val[1]],[wall.points[1].val[2],wall.points[4].val[2]],linecolor = col,xlims=(0,x_max),ylims=(0,y_max),grid=false,legend=false,axis=false)
  #   end
  #   return town_plan
  # end

  # function plot_walls!(project; ids = [], col = :black)
  #   if ids == []
  #       ids = 1:length(project.plan.walls)
  #   end
  #
  #   for wall in project.plan.walls[ids]
  #     if wall.polygon[1][3]==wall.polygon[3][3]
  #       continue
  #     end
  #     xs = [wall.polygon[1][1],wall.polygon[4][1]]
  #     ys = [wall.polygon[1][2],wall.polygon[4][2]]
  #     plot!(xs,ys,
  #         linecolor=col,
  #         xlims = project.plan.limits[1,:],
  #         ylims = project.plan.limits[2,:],
  #         legend = false)
  #   end
  # end

  # function read_data_2d(data_file)
  #   walls = []
  #   input = open(data_file)
  #   for (line_ind,line) in enumerate(eachline(input))
  #     strvec = split(line[1:end-2]," ")
  #     v = map(x->parse(Float64,x),strvec)
  #     wall = reshape(v[:],2,Int(length(v)/2))'
  #     push!(walls,wall)
  #   end
  #   return walls
  # end


  # function create_all_vertex_paths(polygon1::Array{Array{Float64}},polygon2::Array{Array{Float64}})
  #   paths = Array{Line}(length(polygon1)*length(polygon2))
  #   path_counter = 1
  #   for i = 1:length(polygon1)
  #     for j = 1:length(polygon2)
  #       paths[path_counter] = Line(polygon1[i],polygon2[j])
  #       path_counter += 1
  #     end
  #   end
  #   return paths
  # end


  function centerDistance(pol1,pol2)
    # calculate the centers of two polygons and find the distance between them
    c1 = sum(pol1)/length(pol1)
    c2 = sum(pol2)/length(pol2)
    return norm(c1-c2)
  end

  function polygonCenter(pol)
      return sum(pol)/length(pol)
  end


  function shortestPolygonDistance(pol1,pol2)
    # iterate over all possible pairs of polygon vertex and find the minimal
    # distance
    min_dist  = 1.e50
    for i = 1:length(pol1),j = 1:length(pol2)
        min_dist = min(min_dist,norm(pol1[i]-pol2[j]))
    end
    return min_dist
  end


  function apVisibIndex(plan::mapPlan,AP)
    # Creates a vector of bool that indicates that corresponding wall is visible
    # from the standpoint of an AP
    ap_visibility = Array{Bool}(length(plan.walls))*false
    ap_list = push!(Array{Array{Float64}}(0),AP)

    for (couple_ind,couple_wall) in enumerate(plan.walls)
      shrkd_wll_plgn2 = shrink_polygon(couple_wall.polygon,.2)

      if centerDistance(ap_list,shrkd_wll_plgn2) > 200.
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


  function checkPolygonVisibility(pol1,pol2,plan,filter)::Bool
    # Return true if there is a LOS between some vertex of two polygons
    # there is a possibility to make this function concurrent
    for v1 = pol1, v2 = pol2
      if no_walls_on_path(Line(v1,v2),
                          plan,
                          filter)
        return true
      end
    end
    return false
  end

  # function wallVisIndex(plan::mapPlan)::Array{Bool}
  #   # creates a matrix of Bool that indictes visivility of a particular wall
  #   # from the standpoint of another wall
  #   function initfc_bool(s::SharedArray)
  #       for i = eachindex(s)
  #           s[i] = false
  #       end
  #   end
  #   now = length(plan.walls)
  #   vis_matr = SharedArray{Bool}((now,now),init = initfc_bool)
  #   print("Inspecting walls....")
  #   @sync @parallel for wall_id = 1:now
  #   #   print("\rInspecting wall $(wall_id)/$(length(plan.walls))...    ")
  #     if plan.walls[wall_id].special
  #         # if wall is special fill corresponding row with true
  #         for couple_id = 1:now
  #           vis_matr[wall_id,1] = true
  #         end
  #     else
  #         # otherwise check visibility
  #         # Shring walls to prevent false paths in corners.
  #         wall1 = shrink_polygon(plan.walls[wall_id].polygon,0.8)
  #         # this visibility check analyzes corners only. To add visibility for long
  #         # walls add center of the wall to the analysis
  #         push!(wall1,polygonCenter(wall1))
  #         # Filter target walls from the query
  #         filter = Array{Int}(2);
  #         # Test against walls that are above matrix diagonal.
  #         for couple_id = (wall_id+1):now
  #           wall2 = shrink_polygon(plan.walls[couple_id].polygon,0.8)
  #           push!(wall2,polygonCenter(wall2))
  #           # Update filter values here to reduce memory usage
  #           filter[1] = wall_id; filter[2] = couple_id;
  #           result = checkPolygonVisibility(wall1,
  #                                           wall2,
  #                                           plan,
  #                                           filter)
  #           vis_matr[wall_id,couple_id] = result
  #       end
  #     end
  #   end
  #
  #   for i = 1:now
  #     if vis_matr[i,1]
  #         for j = 1:now
  #           vis_matr[i,j] = true
  #           vis_matr[j,i] = true
  #         end
  #     else
  #         for j = (i+1):now
  #           vis_matr[j,i] = vis_matr[i,j]
  #         end
  #     end
  #     # the wall should not be visible from its own standpoint
  #     vis_matr[i,i] = false
  #   end
  #
  #   print("\n")
  #
  #   return vis_matr
  # end



  function wallVisIndexProj(plan::mapPlan)::Array{Bool}
    # creates a matrix of Bool that indictes visivility of a particular wall
    # from the standpoint of another wall
    function initfc_bool(s::SharedArray)
        for i = eachindex(s)
            s[i] = false
        end
    end
    now = length(plan.walls) # Number of Walls
    vis_matr = SharedArray{Bool}((now,now),init = initfc_bool)
    # print("Inspecting walls....")
    for wall_id = 1:now
        print("\rInspecting walls $wall_id....")
        # Get polygon for the first wall
        wall_pol = plan.walls[wall_id].polygon
        xy_p = shrink_line!(Line(wall_pol[1][1:2],wall_pol[4][1:2]),-0.8)
        yz_p = shrink_line!(Line(wall_pol[1][2:3],wall_pol[2][2:3]),-0.8)
        wall1 = shrink_polygon(plan.walls[wall_id].polygon,0.8)
        push!(wall1,polygonCenter(wall1))
        filter = Array{Int}(2);



        @sync @parallel for couple_id = (wall_id+1):now
            # Get polygon for the second wall
            wall_pol2 = plan.walls[couple_id].polygon
            xy_p2 = shrink_line!(Line(wall_pol2[1][1:2],wall_pol2[4][1:2]),-0.8)
            yz_p2 = shrink_line!(Line(wall_pol2[1][2:3],wall_pol2[2][2:3]),-0.8)

            result = false

            # Checking of walls overlap in 2 projections
            # Walls are assumed to be rectangular
            if lines_crossed(xy_p,xy_p2)
                result = true
            elseif lines_crossed(yz_p,yz_p2)
                result = true
            end

            if !result
                wall2 = shrink_polygon(wall_pol2,0.8)
                push!(wall2,polygonCenter(wall2))
                # Update filter values here to reduce memory usage
                filter[1] = wall_id; filter[2] = couple_id;
                result = checkPolygonVisibility(wall1,
                                                wall2,
                                                plan,
                                                filter)
            end
            vis_matr[wall_id,couple_id] = result
        end
    end

    for i = 1:now
      for j = (i+1):now
        vis_matr[j,i] = vis_matr[i,j]
      end
      # the wall should not be visible from its own standpoint
      vis_matr[i,i] = false
    end

    print("\n")

    return convert(Array{Bool},vis_matr)
  end



  function no_walls_on_path(path::Line,plan::mapPlan,filter)::Bool
    shrink_line!(path,float_err_marg)
    for wall_id = query_walls(path,plan.index)
      if MapPrimitives.get_intersection_point(path,plan.walls[wall_id])!=-1
        if wall_id in filter
          continue
        end
        return false
      end
    end
    return true
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
