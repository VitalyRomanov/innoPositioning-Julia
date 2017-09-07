module ImageTree

  using Geometry
  using MapPrimitives
  using MapPlan


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
        # when considering wall visibility from the position of AP
        # we can instantly decide what is visible
      feasible_walls = plan.walls[find(AP_visibility)]
    else
        # otherwise we need to use two step filtering
        # first we find walls that are visible from the posiotion of reflection
        # wall
      visible_walls = plan.walls[find(plan.vis_matr[image.assigned_wall,:])]

    #   then we divide the space into two regions: on both sides of
    #   reflection wall
    #   define the reference position
      parent_wall = image_tree[image.parent].assigned_wall
      if parent_wall == -1
        parent_ref = image_tree[1].location
      else
        parent_ref = point_on_wall(plan.walls[parent_wall])
      end
    #   we say that the target wall should be located on the same side
    #   of the reflection wall as the parent wall
      feasible_walls = walls_on_the_same_side(visible_walls,parent_ref,plan.walls[image.assigned_wall])
      # current_wall = plan.walls[image.assigned_wall]
      # parent_position = normal_directed_towards_point(parent_reference,current_wall)
    end
    return feasible_walls
  end


  function sanityCheck(image_tree,image_id)
    # returns false if there are two reflections from the wall with image_id
    # created to check for multiple reflections from floor (impossible)
    image = image_tree[image_id]
    while image.parent!=-1
      current_wall = image.assigned_wall
      if current_wall == 1 # 1 is the id of floor
        return false
      end
      image = image_tree[image.parent]
    end
    return true
  end


  function calculate_offsprings(image_tree::Array{treeNode},image_id::Int,plan::mapPlan,AP_visibility;max_levels = 3, dist_thr = 150)
    # This function calculates the offsprings of a particular image in the image tree.
    # It is crucial to reduce the number of irrelevant images as much a possible. For this reason several optimizations are done.
    # 1. If the node is the root node, the walls are checked for being visible from the standpoint of the access point.
    # 2. Otherwise visibility matrix is used to determine whether the wall is visible from the standpoint of reflection wall
    # 3. The images should be built only for walls that are on the same side as the previous secondary source
    tree_size = length(image_tree)
    offsprings = Array{treeNode}(length(plan.walls))
    children = Array{Int}(length(plan.walls))
    real_offsprings = 0
    current_offset = 0

    image = image_tree[image_id]

    if image.level < max_levels
      # if the image is the root image no prior knowledge of wall visibility is available
      feasible_walls = get_feasible_walls(image_tree,image,plan,AP_visibility)

      for wall in feasible_walls
        if wall.id == 1
          if !sanityCheck(image_tree,image_id)
            continue
          end
        end
        position = reflection_from_plane(image.location,wall.plane_eq)
        new_node = ImageTree.treeNode(position,image.level+1,image_id,Array{Int}(0),wall.id,norm(position-image.location)/2)
        real_offsprings += 1
        current_offset += 1
        offsprings[real_offsprings] = new_node
        children[real_offsprings] = tree_size+current_offset
      end
    end

    image.children = children[1:real_offsprings]
    return offsprings[1:real_offsprings]
  end



  function build_image_tree(plan::mapPlan,AP::Array{Float64},AP_visibility; pl_thr_dist = 180.)
    image_tree = Array{treeNode}(0)

    root = treeNode(AP,0,-1,[],-1,0.)
    push!(image_tree,root)

    println(plan.vis_matr[:,1])
    println(plan.vis_matr[1,:])

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
  end

end
