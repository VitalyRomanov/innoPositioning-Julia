module ImageTree

  using MapPrimitives


  export treeNode,calculate_offsprings,max_levels
  type treeNode
    location::Point
    level::Int
    parent::Int
    children::Array{Int}
    assigned_wall::Int
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
