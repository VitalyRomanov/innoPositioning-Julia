module Geometry

  const float_err_marg = 0.001

  export reflection_from_plane,get_direction_vector,point_is_on_line,point_is_in_rectangle,shrink_line!,line2mbr,line_plane_intersection,lines_crossed,shrink_polygon,enumerate_mbr2d_geometry
  export Line,MBR,float_err_marg


  type Line
    v1::Array{Float64}
    v2::Array{Float64}
  end

  Line() = Line([0.,0,0],[0.,0,0])

  type MBR
    v1::Array{Float64}
    v2::Array{Float64}
  end

  function line2mbr(line::Line)
    temp = [line.v1 line.v2]
    v1 = minimum(temp,2)
    v2 = maximum(temp,2)
    return MBR(v1,v2)
  end

  # @views (function reflection_from_plane(point,plane)
  #   # make sure the reflection is on the side of the plane that is closer to the second vertex
  #   cp = (dot(point,plane[1:3])+plane[4])/dot(plane[1:3],plane[1:3])
  #   # cp = dot([point;1],plane)/norm(plane[1:3])^2
  #   return point-2*plane[1:3].*cp
  # end)
  @views (function reflection_from_plane(point,plane)
    temp = plane[1:3]
    cp = (dot(point,temp).+plane[4])./norm(temp).^2
    return point.-2.*temp.*cp
  end)

  function shrink_line!(line::Line,asv::Float64)
    eps = line.v2-line.v1
    lnorm = norm(eps)
    if lnorm==0
      return line
    end
    eps = asv*eps/lnorm
    line.v2 -= eps
    line.v1 += eps
    return line
  end

  function shrink_polygon(pol::Array{Array{Float64}},absolute_shrinkage_value)
    polygon = Array{Array{Float64}}(0)
    if length(pol)>0
      center = 0
      for v in pol
        center += v
      end
      center /= length(pol)
      for ind = 1:length(pol)
        eps = center - pol[ind]
        denom = norm(eps)
        if denom==0||absolute_shrinkage_value>denom
          return pol
        end
        eps = absolute_shrinkage_value.*eps./denom
        push!(polygon,pol[ind] + eps)
      end
    end
    return polygon
  end

  # function enumerate_mbr2d_geometry(mbr2d::MBR)
  #   paths = Array(Line,0)
  #   push!(paths,Line(mbr2d.v1,[mbr2d.v1[1],mbr2d.v2[2]]))
  #   push!(paths,Line([mbr2d.v1[1],mbr2d.v2[2]],mbr2d.v2))
  #   push!(paths,Line(mbr2d.v2,[mbr2d.v2[1],mbr2d.v1[2]]))
  #   push!(paths,Line([mbr2d.v2[1],mbr2d.v1[2]],mbr2d.v1))
  #   return paths
  # end



  function get_direction_vector(v1::Array{Float64},v2::Array{Float64})::Array{Float64}
    return v2-v1
  end

  @inline function get_direction_vector(line::Line)::Array{Float64}
    return line.v2-line.v1
  end


  function line_plane_intersection(line::Line,plane::Array{Float64})
    direction_vector = Geometry.get_direction_vector(line)
    t2 = dot(direction_vector,plane[1:3])
    if t2==0
      return -1
    end
    t1 = dot(line.v1,plane[1:3])
    c = (-plane[4]-t1)/t2
    if abs(c) > Geometry.float_err_marg
      return line.v1+direction_vector.*c
    else
      return -1
    end
  end


  function point_is_on_line(line::Line,ip::Array{Float64})
    return dot(line.v1-ip,line.v2-ip)<=0
  end


  function point_is_in_rectangle(ip::Array{Float64}, polygon::Array{Array{Float64}})
    dv12 = get_direction_vector(polygon[1],polygon[2])
    dv14  = get_direction_vector(polygon[1],polygon[4])
    if dot(polygon[1],dv12) <=
      dot(ip,dv12) <=
      dot(polygon[2],dv12) &&
      dot(polygon[1],dv14) <=
      dot(ip,dv14) <=
      dot(polygon[4],dv14)
      return true
    end
    return false
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

    return dot(line1.v1-ip,line1.v2-ip)<=0 && dot(line2.v1-ip,line2.v2-ip)<=0
  end

end
