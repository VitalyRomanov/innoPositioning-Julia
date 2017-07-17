module MapPrimitives

  using Geometry

  export Point, Wall3D, Sector, Wall
  export get_plane_equation!,wall2mbr,get_intersection_point

  type Wall3D
    id::Int
    polygon::Array{Array{Float64}}
    plane_eq::Array{Float64}
  end

  type Point
    val::Array{Float64}
  end

  type Wall
    id::Int
    points::Array{Point}
    plane_equation::Array{Float64}
  end

  type Sector
    id::Int
    walls::Array{Int}
    geometry::Array{Line}
  end


  # function wall2mbr(wall::Wall)
  #   v1 = minimum([wall.points[1].val wall.points[2].val wall.points[3].val wall.points[4].val],2)
  #   v2 = maximum([wall.points[1].val wall.points[2].val wall.points[3].val wall.points[4].val],2)
  #   return MBR(v1,v2)
  # end

  function wall2mbr(wall::Wall3D)
    v1 = minimum([wall.polygon[1] wall.polygon[2] wall.polygon[3] wall.polygon[4]],2)
    v2 = maximum([wall.polygon[1] wall.polygon[2] wall.polygon[3] wall.polygon[4]],2)
    return MBR(v1,v2)
  end

  # type WallIndex
  #   walls::Array{Wall}
  #   sectors::Array{Sector}
  #   sector_size::Float64
  #   sectors_dim::Array{Int}
  #   space_size::Array{Int}
  # end

 # function get_sector_ind(coord,sector_size,sectors)
 #   x = Int(ceil(coord[1]/sector_size))
 #   y = Int(ceil(coord[2]/sector_size))
 #   if x<1 || y<1 || x >= sectors[1] || y >= sectors[2]
 #     return -1
 #   end
 #  #  print(coord[1]," ",coord[2]," ",x," ",y," ",sectors[1]," ",sectors[2],"\n")
 #   return (y-1)*sectors[1]+x
 # end
 #
 # function get_2d_sector_ind(coord,sector_size,sectors)
 #   x = Int(ceil(coord[1]/sector_size))
 #   y = Int(ceil(coord[2]/sector_size))
 #   return [x,y]
 # end




  # function create_wall_index(walls::Array{Wall},x_max,y_max,sector_size = 30.0)
  #   x_max_grid = Int64(ceil(x_max/sector_size))
  #   y_max_grid = Int64(ceil(y_max/sector_size))
  #   number_of_sector = x_max_grid*y_max_grid
  #   wall_index = WallIndex(walls,Array(Sector,number_of_sector),sector_size,[x_max_grid,y_max_grid],[x_max,y_max])
  #   # print(prod(wall_index.sectors_dim),"\n")
  #
  #   sector_count = 1
  #   for y_grid in 1:y_max_grid
  #     for x_grid in 1:x_max_grid
  #       new_sector = Sector(sector_count,Array(Int,0),Array(Line,0))
  #       push!(new_sector.geometry,Line([(x_grid-1)*sector_size,(y_grid-1)*sector_size],[(x_grid)*sector_size,(y_grid-1)*sector_size]))
  #       push!(new_sector.geometry,Line([(x_grid)*sector_size,(y_grid-1)*sector_size],[(x_grid)*sector_size,(y_grid)*sector_size]))
  #       push!(new_sector.geometry,Line([(x_grid)*sector_size,(y_grid)*sector_size],[(x_grid-1)*sector_size,(y_grid)*sector_size]))
  #       push!(new_sector.geometry,Line([(x_grid-1)*sector_size,(y_grid)*sector_size],[(x_grid-1)*sector_size,(y_grid-1)*sector_size]))
  #       wall_index.sectors[sector_count] = new_sector
  #       sector_count += 1
  #     end
  #   end
  #
  #   # for i in 1:prod(wall_index.sectors_dim)
  #   #   push!(wall_index.sectors,Sector(i,Array(Int,0),Array(Wall,0)))
  #   #   push!(wall_index.sectors[i].walls,length(walls))
  #   #   new_wall = Line(,)
  #   #   push!(wall_index.sectors[i].geometry,)
  #   # end
  #
  #   for wall in walls
  #     cs1 = get_sector_ind(wall.points[1].val,wall_index.sector_size,wall_index.sectors_dim)
  #     cs2 = get_sector_ind(wall.points[4].val,wall_index.sector_size,wall_index.sectors_dim)
  #     if cs1==-1 || cs2==-1
  #       continue
  #     end
  #     # print(cs1," ",cs2,"\n")
  #     if cs1==cs2
  #       push!(wall_index.sectors[cs1].walls,wall.id)
  #     else
  #       push!(wall_index.sectors[cs1].walls,wall.id)
  #       push!(wall_index.sectors[cs2].walls,wall.id)
  #     end
  #   end
  #
  #   # print(wall_index,"\n")
  #   # for wall in walls
  #
  #   return wall_index
  #
  # end


  # type SpaceInfo
  #   space_size
  #   grid
  #   grid_size::Float64
  # end
  #
  # function getSpaceInfo(space_size,grid_size=1.)
  #   grid = convert(Array{Int32,2},ceil(space_size/grid_size))
  #   return SpaceInfo(space_size,grid,grid_size)
  # end


  function get_plane_equation(cw::Wall)
    equation = [.0,.0,.0,.0]
    equation[1] = (cw.points[2].val[2]-cw.points[1].val[2])*(cw.points[3].val[3]-cw.points[1].val[3])-(cw.points[3].val[2]-cw.points[1].val[2])*(cw.points[2].val[3]-cw.points[1].val[3])
    equation[2] = -((cw.points[2].val[1]-cw.points[1].val[1])*(cw.points[3].val[3]-cw.points[1].val[3])-(cw.points[3].val[1]-cw.points[1].val[1])*(cw.points[2].val[3]-cw.points[1].val[3]))
    equation[3] = (cw.points[2].val[1]-cw.points[1].val[1])*(cw.points[3].val[2]-cw.points[1].val[2])-(cw.points[3].val[1]-cw.points[1].val[1])*(cw.points[2].val[2]-cw.points[1].val[2])
    equation[4] = -cw.points[1].val[1]*equation[1]-cw.points[1].val[2]*equation[2]-cw.points[1].val[3]*equation[3]
    cw.plane_equation = equation
    return equation
  end

  function get_plane_equation!(cw::Wall3D)
    equation = [.0,.0,.0,.0]
    equation[1] = (cw.polygon[2][2]-cw.polygon[1][2])*(cw.polygon[3][3]-cw.polygon[1][3])-(cw.polygon[3][2]-cw.polygon[1][2])*(cw.polygon[2][3]-cw.polygon[1][3])
    equation[2] = -((cw.polygon[2][1]-cw.polygon[1][1])*(cw.polygon[3][3]-cw.polygon[1][3])-(cw.polygon[3][1]-cw.polygon[1][1])*(cw.polygon[2][3]-cw.polygon[1][3]))
    equation[3] = (cw.polygon[2][1]-cw.polygon[1][1])*(cw.polygon[3][2]-cw.polygon[1][2])-(cw.polygon[3][1]-cw.polygon[1][1])*(cw.polygon[2][2]-cw.polygon[1][2])
    equation[4] = -cw.polygon[1][1]*equation[1]-cw.polygon[1][2]*equation[2]-cw.polygon[1][3]*equation[3]
    cw.plane_eq = equation
  end


  function get_intersection_point(line::Line,wall::Wall3D)
    direction_vector = get_direction_vector(line)
    ip = line_plane_intersection(line,wall.plane_eq)
    if ip!=-1
    # cp = plane.plane_eq[1:3]
    # num = dot(cp,(plane.polygon[1]-line.v1))
    # denom = dot(cp,direction_vector)
    # if denom == 0.
    #   return -1
    # end
    # t = num / denom
    # if t > float_err_marg
    #   ip = line.v1+direction_vector.*t
      if point_is_on_line(line,ip)
        if point_is_in_rectangle(ip, wall.polygon)
          return ip
        end
      end
    end
    return -1
  end

end
