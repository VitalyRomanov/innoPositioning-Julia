module MapBuilder
  using MapPrimitives
  using MapPlan
  using ImageTree
  using Geometry
  using RadixTree


  export build_image_tree,calculate_signal_strength,get_visibility_matrix

  function get_path_info(image_tree::Array{treeNode},
                        node_ind::Int,
                        Rx::Array{Float64},
                        plan::mapPlan)
    image = image_tree[node_ind]

    last_vertex = Rx
    wop = 0::Int
    nrw = image.level
    dist = 0.::Float64
    susp_path = Line()

    while image.parent!=-1
      current_image_loc = image.location
      current_wall = image.assigned_wall

      susp_path.v1 = last_vertex; susp_path.v2 = current_image_loc
    #   susp_path = Line(last_vertex,current_image_loc)

      reflection_point = get_intersection_point(susp_path,plan.walls[current_wall])

      if reflection_point != -1
        susp_path.v2 = reflection_point
        # susp_path = Line(last_vertex,reflection_point)
        wop += walls_on_path(susp_path,plan)
        dist += norm(last_vertex - reflection_point)
      else
        return -1.,nrw,wop
      end

      last_vertex = reflection_point
      image = image_tree[image.parent]
    end

    susp_path.v1 = last_vertex; susp_path.v2 = image.location
    # susp_path = Line(last_vertex,image.location)
    wop += walls_on_path(susp_path,plan)
    dist += norm(last_vertex - image.location)

    return dist,nrw,wop
  end



  function get_path_vertices(image_tree::Array{treeNode},
                            node_ind::Int,
                            Rx::Array{Float64},
                            plan::mapPlan)

    image = image_tree[node_ind]
    vertices = Array{Array{Float64}}(image.level+2)
    vertices[end] = Rx
    vert_pos = length(vertices)-1

    while image.parent!=-1

      susp_path = Line(vertices[vert_pos+1],image.location)
      intersection_point = get_intersection_point(susp_path,plan.walls[image.assigned_wall])

      if intersection_point == -1
        return Array{Array{Float64}}(0)
      end

      vertices[vert_pos] = intersection_point
      vert_pos -= 1
      image = image_tree[image.parent]
    end

    vertices[vert_pos] = image.location

    return vertices
  end

  function dump!(fileobj,d,r,i,x,y)
    for ind = 1:length(d)
      if d[ind]!=-1.
        write(fileobj,"$x $y $(d[ind]) $(r[ind]) $(i[ind])\n")
      end
    end
  end

  function getDumpPos(filePath)
    dumpfile = open(filePath,"r")
    lines = readlines(dumpfile)
    if length(lines)>0
      x,y = map(x->parse(Float64,x),split(lines[end]," "))[1:2]
    else
      x,y = (1,1)
    end
    return Int(x+1),Int(y+1)
  end

  function caclulate_signal_strength_matrix(image_tree::Array{treeNode},
                                            plan::mapPlan,
                                            parameters::Array{Float64},
                                            dumpPath)
    # add support for multiple resolutions
    size1 = plan.limits[1,2]-plan.limits[1,1]; size2 = plan.limits[2,2]-plan.limits[2,1];
    x_start = 1; y_start = 1;

    if isfile(dumpPath)
      x_start,y_start = getDumpPos(dumpPath)
    end

    dumpfile = open(dumpPath,"a+")

    # ss = SharedArray{Float64}(size1,size2)
    ss = Array{Float64}(size1,size2)
    for x = x_start:size1
      # @sync @parallel for y = 1:size2
      for y = y_start:size2
        # incorporate cell size of type Float64
        x_loc = Float64(x + plan.limits[1,1]) # multiply by cell size
        y_loc = Float64(y + plan.limits[2,1]) # multiply by cell size
        ss[x,y],d,r,i = calculate_signal_strength([x_loc,y_loc,1],image_tree,plan,parameters)
        dump!(dumpfile,d,r,i,x,y)
        print("\r$(x)/$(size1) $(y)/$(size2) $(ss[x,y])      ")
      end
      # print("\r$(x)/$(size1)       ")
      y_start = 1
    end
    print("\n")
    close(dumpfile)
    return ss
  end



  function signal_paths_info(Rx::Array{Float64},image_tree::Array{treeNode},plan::mapPlan; ae = 2,rc = 12.53,tc = 5., P0 = 0., sc = 0.5, pldthr = 200.)
    # for each path three parameters are estimated: 10*log10(path_dist),
    # number_of_reflections, number_of_walls_on_path
    paths_info = Array{Float64}(0,3)

    for (node_ind,image) in enumerate(image_tree)
      if norm(image.location - Rx) > pldthr
        continue
      end

      dist,nrw,wop = get_path_info(image_tree,node_ind,Rx,plan)

      if dist!=-1.
        paths_info = [paths_info;[10*log10(dist) nrw wop]]
      end

    end

    return paths_info
  end


  function calculate_signal_strength(Rx::Array{Float64},
                                    image_tree::Array{treeNode},
                                    plan::mapPlan,
                                    parameters::Array{Float64};
                                    pldthr = 150.)

    # This function initializes a SharedArray to zero
    function initfn(S::SharedArray)
      for i in eachindex(S)
        S[i] = 0.
      end
    end

    # Shared array should be initialized in a "SharedArray way"
    distan = SharedArray{Float64}(length(image_tree),1, init = initfn)
    reflec = SharedArray{Float64}(length(image_tree),1, init = initfn)
    inters = SharedArray{Float64}(length(image_tree),1, init = initfn)
    path_loss = SharedArray{Float64}(length(image_tree),1, init = initfn)

    used_images = 0

    equation_line = [1. 1 1 1 0 0 0]

    @sync @parallel for image_id = 1:length(image_tree)
      # !!!!!!!!!!!!!! the following conditions is suspected to be wrong
      image = image_tree[image_id]
      if norm(image.location - Rx) < pldthr

        dist,nrw,wop = get_path_info(image_tree,image_id,Rx,plan)
        distan[image_id] = dist
        reflec[image_id] = nrw
        inters[image_id] = wop

        #   this is intended to skip paths that go through walls
        if wop==0 && dist!=-1.
          used_images += 1
          equation_line[5] = log10(dist)
          equation_line[6] = nrw
          equation_line[7] = wop

          pl = 10 ^ ( sum(equation_line*parameters) / 10. )
          # pl = 10 ^ ( sum([1 1 1 1 log10(dist) nrw wop]*parameters) / 10. )
          # pl = 10^( ( P0 - (10*ae*log10(dist)+20*log10(2.4e9) - 147.55 + sc*rc*nrw + tc*wop)) /10.)
          # push!(path_loss,pl)
          path_loss[image_id] = pl
        end
      end
    end

    pl_sum = sum(path_loss)
    if pl_sum==0
      value = -900.
    else
      value = 10*log10(pl_sum)
    end

    return value,distan,reflec,inters
  end




  function calculate_paths(Rx::Array{Float64},
                          image_tree::Array{treeNode},
                          plan::mapPlan;
                          pldthr = 200.)

    paths = []

    for (node_ind,image) in enumerate(image_tree)
      # !!!!!!!!!!!!!! the following conditions is suspected to be wrong
      if norm(image.location - Rx) > pldthr
        continue
      end

      vert = get_path_vertices(image_tree,node_ind,Rx,plan)

      if length(vert)>0
        push!(paths,vert)
      end

    end

    return paths
  end

  # function get_path_distance(path)
  #   dist = 0.
  #   for i in 1:(length(path)-1)
  #     dist += norm(path[i]-path[i+1])
  #   end
  #   return dist
  # end


end
