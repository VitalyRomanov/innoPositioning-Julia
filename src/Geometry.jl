module Geometry

  using MapPrimitives

  export reflection_from_plane,get_direction_vector,get_intersection_point,check_if_point_belongs_to_line_segment,check_if_point_belongs_to_wall

  function reflection_from_plane(point,plane)
    cp = dot([point;1],plane)/norm(plane[1:3])^2
    return point-2*plane[1:3].*cp
  end

  function get_direction_vector(v1,v2)
    return v2-v1
  end

  function get_intersection_point(v1,v2,plane)
    intersect_epsilon = 0.1
    direction_vector = get_direction_vector(v1,v2)
    cp = plane.plane_equation[1:3]
    num = dot(cp,(plane.points[1].val-v1))
    denom = dot(cp,direction_vector)
    if denom == 0.
      return -1
    end
    t = num / denom
    if t > intersect_epsilon
      ip = v1+direction_vector.*t
      if check_if_point_belongs_to_line_segment(v1,v2,ip)
        if check_if_point_belongs_to_wall(ip, plane)
          return ip
        end
      end
    end
    return -1
  end

  function check_if_point_belongs_to_line_segment(v1,v2,ip)
    return dot(v1-ip,v2-ip)<=0
  end


  function check_if_point_belongs_to_wall(ip, plane)
    dv12 = get_direction_vector(plane.points[1].val,plane.points[2].val)
    dv14  = get_direction_vector(plane.points[1].val,plane.points[4].val)
    if dot(plane.points[1].val,dv12) <=
      dot(ip,dv12) <=
      dot(plane.points[2].val,dv12) &&
      dot(plane.points[1].val,dv14) <=
      dot(ip,dv14) <=
      dot(plane.points[4].val,dv14)
      return true
    end
    return false
  end

end
