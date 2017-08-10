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


  function wall2mbr(wall::Wall3D)
    v1 = minimum([wall.polygon[1] wall.polygon[2] wall.polygon[3] wall.polygon[4]],2)
    v2 = maximum([wall.polygon[1] wall.polygon[2] wall.polygon[3] wall.polygon[4]],2)
    return MBR(v1,v2)
  end


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
      if point_is_on_line(line,ip)
        if point_is_in_rectangle(ip, wall.polygon)
          return ip
        end
      end
    end
    return -1
  end

end
