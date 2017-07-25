module ImageTree

  using Geometry
  using MapPrimitives
  using MapPlan
  # using RadixTree


  export treeNode,calculate_offsprings,build_image_tree
  type treeNode
    location::Array{Float64}
    level::Int
    parent::Int           # equals -1 for the root node
    children::Array{Int}
    assigned_wall::Int
    lb_dist::Float64
  end


  function walls_on_the_same_side(walls, ref_wall_point, refl_wall)
    # plane_eq = refl_wall.plane_eq
    temp = normal_directed_towards_point(ref_wall_point,refl_wall)
    correct_side = map(x->(temp==normal_directed_towards_point(wall_center(x),refl_wall)),walls)
    return walls[find(correct_side)]
  end

  function wall_center(wall::Wall3D)
    center = 0
    for v in wall.polygon
      center += v
    end
    center /= length(wall.polygon)
    return center
  end

  function point_on_wall(wall::Wall3D)
    return wall.polygon[1]
  end


  function get_feasible_walls(image_tree,image,plan,AP_visibility)
    if image.assigned_wall == -1
      feasible_walls = plan.walls[find(AP_visibility)]
    else
      visible_walls = plan.walls[find(plan.vis_matr[image.assigned_wall,:])]

      parent_wall = image_tree[image.parent].assigned_wall
      if parent_wall == -1
        parent_ref = image_tree[1].location
      else
        parent_ref = point_on_wall(plan.walls[parent_wall])
      end
      feasible_walls = walls_on_the_same_side(visible_walls,parent_ref,plan.walls[image.assigned_wall])
      # current_wall = plan.walls[image.assigned_wall]
      # parent_position = normal_directed_towards_point(parent_reference,current_wall)
    end
    return feasible_walls
  end

  function calculate_offsprings(image_tree::Array{treeNode},image_id::Int,plan::mapPlan,AP_visibility::Array{Bool};max_levels = 3, dist_thr = 150)
    # This function calculates the offsprings of a particular image in the image tree.
    # It is crucial to reduce the number of irrelevant images as much a possible. For this reason several optimizations are done.
    # 1. If the node is the root node, the walls are checked for being visible from the standpoint of the access point.
    # 2. Otherwise visibility matrix is used to determine whether the wall is visible from the standpoint of reflection wall
    # 3. The images should be built only for walls that are on the same side as the previous secondary source
    tree_size = length(image_tree)
    offsprings = Array(treeNode,length(plan.walls))
    children = Array(Int,length(plan.walls))
    real_offsprings = 0
    current_offset = 0

    image = image_tree[image_id]

    if image.level < max_levels
      # if the image is the root image no prior knowledge of wall visibility is available
      feasible_walls = get_feasible_walls(image_tree,image,plan,AP_visibility)
      # if image.assigned_wall == -1
      #   feasible_walls = plan.walls[find(AP_visibility)]
      # else
      #   visible_walls = plan.walls[find(plan.vis_matr[image.assigned_wall,:])]
      #
      #   parent_wall = image_tree[image.parent].assigned_wall
      #   if parent_wall == -1
      #     parent_ref = image_tree[1].location
      #   else
      #     parent_ref = point_on_wall(plan.walls[parent_wall])
      #   end
      #   feasible_walls = walls_on_the_same_side(visible_walls,parent_ref,plan.walls[image.assigned_wall])
      #   # current_wall = plan.walls[image.assigned_wall]
      #   # parent_position = normal_directed_towards_point(parent_reference,current_wall)
      # end

      # we are going to reduce the number of feasible walls
      # temp_feasible_walls = Array(Wall3D,length(feasible_walls))
      # temp_position = 1

      # the distance between the original image and the current image is not a very good measure of feasibility
      # instead we can calculate the sum distance between consecutive images that will give the distance lowe bound
      # for wall in feasible_walls
      #   if norm(wall.polygon[1]-image_tree[1].location)<distance_threshold
      #     if image.assigned_wall == -1
      #       truly_feasible =  no_walls_on_path(Line(image.location,(wall.polygon[1]+wall.polygon[4])/2),plan)
      #     else
      #       truly_feasible = (parent_position==normal_directed_towards_point(wall.polygon[1],current_wall))
      #     end
      #
      #     if truly_feasible
      #       temp_feasible_walls[temp_position] = wall
      #       temp_position += 1
      #     end
      #   end
      # end

      # feasible_walls = temp_feasible_walls[1:temp_position-1]

      for wall in feasible_walls
        position = reflection_from_plane(image.location,wall.plane_eq)
        new_node = ImageTree.treeNode(position,image.level+1,image_id,Array(Int,0),wall.id,norm(position-image.location)/2)
        real_offsprings += 1
        current_offset += 1
        offsprings[real_offsprings] = new_node
        children[real_offsprings] = tree_size+current_offset
      end
    end

    image.children = children[1:real_offsprings]
    return offsprings[1:real_offsprings]
  end



  function build_image_tree(plan::mapPlan,AP::Array{Float64},AP_visibility::Array{Bool}; pl_thr_dist = 180.)
    image_tree = Array(treeNode,0)

    root = treeNode(AP,0,-1,[],-1,0.)
    push!(image_tree,root)

    tree_position = 1
    iteration_steps = length(image_tree)

    while tree_position <= iteration_steps
      if image_tree[tree_position].lb_dist < pl_thr_dist
        offsprings = calculate_offsprings(image_tree,tree_position,plan,AP_visibility,dist_thr = pl_thr_dist)
        append!(image_tree,offsprings)
        iteration_steps = length(image_tree)
      end
      tree_position += 1
    end

    print("Image tree is built: ",length(image_tree),"\n")
    return image_tree
  end


  function normal_directed_towards_point(test_pos::Array{Float64},wall::Wall3D)
    # this function heeds to be adopted for 3D
    direction_vect = wall.plane_eq[1:3]
    t2 = dot(direction_vect,direction_vect)
    t1 = dot(test_pos,direction_vect)
    c = (-wall.plane_eq[4]-t1)/t2
    pow = test_pos+direction_vect.*c
    return dot(pow-test_pos,direction_vect)<=0
    # slope = (reference_wall.polygon[4][2]-reference_wall.polygon[1][2])/(reference_wall.polygon[4][1]-reference_wall.polygon[1][1])
    # bias = reference_wall.polygon[4][2] - slope * reference_wall.polygon[4][1]
    # return test_pos[2] > slope * test_pos[2] + bias
  end

  # function calculate_offsprings(image::treeNode,image_id::Int,tree_size::Int,walls::Array{Wall},AP::Point,visibility_matrix;max_levels = 3, distance_threshold = 150)
  #   offsprings = Array(treeNode,length(walls))
  #   children = Array(Int,length(walls))
  #   real_offsprings = 0
  #   current_offset = 0
  #
  #   if image.level < max_levels
  #     if image.assigned_wall == -1
  #       feasible_walls = walls
  #     else
  #       feasible_walls = walls[find(visibility_matrix[image.assigned_wall,:])]
  #     end
  #     for wall in feasible_walls
  #       if norm(wall.points[1].val-AP.val)<distance_threshold
  #         if wall.id != image.assigned_wall
  #           paths = Array(Line,2)
  #           paths[1] = Line(image.location.val[1:2],wall.points[1].val[1:2])
  #           paths[2] = Line(image.location.val[1:2],wall.points[4].val[1:2])
  #
  #           wall_visible = false
  #           for path in paths
  #             if no_walls_on_path(path,wall_ind)
  #               wall_visible = true
  #               break
  #             end
  #           end
  #           if wall_visible
  #             cp = dot([image.location.val;1],wall.plane_equation)/norm(wall.plane_equation[1:3])^2
  #             position = image.location.val-2*wall.plane_equation[1:3].*cp
  #             new_image = MapPrimitives.Point(position)
  #             new_node = ImageTree.treeNode(new_image,image.level+1,image_id,Array(Int,0),wall.id)
  #             # push!(image.children,tree_size+current_offset)
  #             # push!(offsprings,new_node)
  #             real_offsprings += 1
  #             current_offset+=1
  #             offsprings[real_offsprings] = new_node
  #             children[real_offsprings] = tree_size+current_offset
  #           end
  #         end
  #       end
  #     end
  #   end
  #   # potentially there can be a problem of not all valuable nodes included in the tree
  #   image.children = children[1:real_offsprings]
  #   return offsprings[1:real_offsprings]
  # end
end
