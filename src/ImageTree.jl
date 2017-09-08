module ImageTree

  using Geometry
  using MapPrimitives
  using MapPlan


  export treeNode,calcu_offsprings,buildImageTree
  type treeNode
    location::Array{Float64}
    level::Int
    parent::Int           # equals -1 for the root node
    children::Array{Int}
    assigned_wall::Int
    lb_dist::Float64
  end


  function walls_on_the_same_side(walls, ref_wall_point, refl_wall)
    # returns walls that are that are on the same side from refleciton wall as
    # the parent wall
    temp = normal_directed_towards_point(ref_wall_point,refl_wall)
    correct_side = map(x->(temp==normal_directed_towards_point(getWallCenter(x),refl_wall)),walls)
    return walls[find(correct_side)]
  end

  function getWallCenter(wall::Wall3D)
    center = 0
    for v in wall.polygon
      center += v
    end
    center /= length(wall.polygon)
    return center
  end

  function getPrevWallCenter(image_tree,image_id,plan)
    # Get the center of a previous a parent wall. If previous wall is special,
    # go one level deeper. If the parent is the root image, return it's location
    center = image_tree[image_id].location

    dive = true
    while dive
        wall_id = image_tree[image_id].assigned_wall
        if wall_id != -1
            if !plan.walls[wall_id].special
                center = getWallCenter(plan.walls[wall_id])
                dive = false
            end
            image_id = image_tree[image_id].parent
        else
            center = image_tree[image_id].location
            dive = false
        end
    end
    return center
  end

  function point_on_wall(wall::Wall3D)
    return wall.polygon[1]
  end


  function getFeasibleWalls(image_tree,image,plan,apVis)
    if image.assigned_wall == -1
        # when considering wall visibility from the position of AP
        # we can instantly decide what is visible
      feasible_walls = plan.walls[find(apVis)]
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


  function sanityCheck(image_tree,image_id,plan)
    # returns false if there are two reflections from special walls
    # created to check for multiple reflections from floor or ceiling
    image = image_tree[image_id]
    while image.parent!=-1
      if plan.walls[image.assigned_wall].special
        return false
      end
      image = image_tree[image.parent]
    end
    return true
  end


  function getOffsprings(image_tree::Array{treeNode},   # tree
                        image_id::Int,                  # current parent node id
                        plan::mapPlan,                  # plan
                        apVis;                          # visibility vector
                        max_levels = 3,                 # do not create new levels beyond that
                        dist_thr = 40.)                 # max dist between walls with conceq reflections
    # This function calculates the offsprings of a particular image in the
    # image tree.
    # It is crucial to reduce the number of irrelevant images as much a
    # possible. For this reason several optimizations are done.
    # 1. If the node is the root node, the walls are checked for being visible
    #   from the standpoint of the access point.
    # 2. Otherwise visibility matrix is used to determine whether the wall is
    #   visible from the standpoint of reflection wall
    # 3. The images should be built only for walls that are on the same side
    #   as the previous secondary source
    tree_size = length(image_tree)
    offsprings = Array{treeNode}(length(plan.walls))
    children = Array{Int}(length(plan.walls))
    real_offsprings = 0
    current_offset = 0

    image = image_tree[image_id]

    if image.level < max_levels
      # if the image is the root image no prior knowledge of wall visibility
      # is available
      feasible_walls = getFeasibleWalls(image_tree,image,plan,apVis)

      for wall in feasible_walls
        # get distance
        if wall.width<1.
            continue
        end
        if image.assigned_wall != -1
            # distBetwWalls = norm(getWallCenter(wall)-getWallCenter(plan.walls[image.assigned_wall]))
            distBetwWalls = norm(getWallCenter(wall)-getPrevWallCenter(image_tree,image_id,plan))
            if (!wall.special)&&(distBetwWalls>dist_thr/image.level)
                continue
            end
        end
        if (wall.special)&&(!sanityCheck(image_tree,image_id,plan))
            continue
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



  function buildImageTree(plan::mapPlan,AP::Array{Float64},apVis; pl_thr_dist = 100.)
    # this function creates an image tree
    # plan:
    # AP: location of the current AP
    # apVis: elements of this vector indicate whether walls are visible form the
    # standpoint of this AP
    # pl_thr_dist: distance between images after which the signal if too week
    # Array for storing tree node
    image_tree = Array{treeNode}(0)

    root = treeNode(AP,0,-1,[],-1,0.)
    push!(image_tree,root)

    # println(plan.vis_matr[:,1])
    # println(plan.vis_matr[1,:])

    # Prepare to populate the tree
    tree_position = 1; iteration_steps = length(image_tree);
    # While new nodes can be inserted, continue
    while tree_position <= iteration_steps
      # Stop if the lower bound distance between images is too large
      if image_tree[tree_position].lb_dist < pl_thr_dist
        offsprings = getOffsprings(image_tree,
                                        tree_position,
                                        plan,
                                        apVis)
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
