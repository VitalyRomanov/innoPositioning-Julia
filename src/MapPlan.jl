module MapPlan
  using MapPrimitives
  using Plots
  export read_data,walls,bbox,AP,wall_ind,plot_walls,downsample_maps

  # function nearby_walls(AP::Point)
  #   loc = get_2d_sector_ind(AP.val,wall_ind.sector_size,wall_ind.sectors_dim)
  #   xmin = loc[1]-Int(ceil(distance_threshold/wall_ind.sector_size))
  #   xmax = loc[1]+Int(ceil(distance_threshold/wall_ind.sector_size))
  #   ymin = loc[2]-Int(ceil(distance_threshold/wall_ind.sector_size))
  #   ymax = loc[2]+Int(ceil(distance_threshold/wall_ind.sector_size))
  #   walls = Set()
  #   for x=xmin:xmax
  #     for y = ymin:ymax
  #       ind = (y-1)*wall_ind.sectors_dim[1]+x
  #       union!(walls,wall_ind.sectors[ind].walls)
  #     end
  #   end
  #   return walls
  # end

  function downsample_maps(signal_maps,factor)
    println("Downsampling maps")
    new_maps = Array(Array{Float64},length(signal_maps))
    for (ind,map) in enumerate(signal_maps)
      new_maps[ind] = downsample_map(map,2)
    end
    return new_maps
  end

  function downsample_map(signal_map,factor)
    new_size = [Int(floor(size(signal_map,1)/factor)),Int(floor(size(signal_map,2)/factor))]
    new_map = Array(Float64,new_size[1],new_size[2])
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


  function read_data(data_file)
    walls=Array(Wall,0)
    bbox=Array(Wall,0)
    input = open(data_file)
    wall_counter = 1
    for (line_ind,line) in enumerate(eachline(input))
      strvec = split(line,",")
      v = map(x->parse(Float64,x),strvec)
      new_wall = Wall(wall_counter,[Point([v[1],v[2],v[3]]),Point([v[4],v[5],v[6]]),Point([v[7],v[8],v[9]]),Point([v[10],v[11],v[12]])],[0,0,0,0])
      MapPrimitives.get_plane_equation(new_wall)
      if line_ind > 4
        wall_counter += 1
        push!(walls,new_wall)
      else
        push!(bbox,new_wall)
      end
    end

    new_wall = Wall(wall_counter+1,[bbox[1].points[1],bbox[2].points[1],bbox[3].points[1],bbox[4].points[1]],[0.,0.,0.,0.])
    MapPrimitives.get_plane_equation(new_wall)
    push!(walls,new_wall)

    print("$(length(walls)) walls imported\n")
    return walls,bbox
  end
end
