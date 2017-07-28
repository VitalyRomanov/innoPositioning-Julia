module RadixTree

using Geometry

const nod = 3 # number of dimensions

# export MBR
#
# type MBR
#   v1::Array{Float64}
#   v2::Array{Float64}
# end



type radixTree
  indexDatasetSize::Int
  numberOfInternalNodes::Int
  objects::Array{MBR}

  id::Array{Int}
  chld::Array{Array{Int}}
  coverage::Array{Array{Int}}
  parent::Array{Int}
  mbr::Array{MBR}
end






function obj2mbr(walls,obj2mbr)
  mbrs = Array(MBR,length(walls))
  for (wall_ind,wall) in enumerate(walls)
    mbrs[wall_ind] = obj2mbr(wall)
  end
  return mbrs
end


function expBits(v::UInt32)
  v = (v * 0x00010001) & 0xFF0000FF
  v = (v * 0x00000101) & 0x0F00F00F
  v = (v * 0x00000011) & 0xC30C30C3
  v = (v * 0x00000005) & 0x49249249
  return v
end

function mortonIn3D(x)
  x = minimum([maximum([x*1024. zeros(Float64,3)],2) ones(Float64,3)*1023.],2)
  xx = expBits(convert(UInt32,round(x[1])))
  yy = expBits(convert(UInt32,round(x[2])))
  zz = expBits(convert(UInt32,round(x[3])))
  return xx * 4 + yy * 2 + zz
end

function get_morton_codes(objects,lims)
  morton_codes = zeros(UInt32,length(objects))

  for (obj_ind,obj) in enumerate(objects)
    coord = ( ( obj.v1 + obj.v2 ) / 2 - lims[:,1] ) ./ lims[:,2]
    morton_codes[obj_ind] = mortonIn3D(coord)
  end

  return morton_codes
end


function determine_range(mortonCodes, i)
  dataSize = length(mortonCodes)
	current = mortonCodes[i]

	dir = leading_zeros(mortonCodes[i]$mortonCodes[i+1])-leading_zeros(mortonCodes[i]$mortonCodes[i-1])
	d = (dir<0)?-1:1

  if i<=6
    println(leading_zeros(mortonCodes[i])," ",mortonCodes[i])
  end

	deltaMin = leading_zeros(current$mortonCodes[i-d])

  lMax = 1
	nextPosition = i+lMax*d
	deltaCurrent = leading_zeros(current$mortonCodes[i+lMax*d])
	while(deltaCurrent > deltaMin)
		lMax *= 2
		nextPosition =  i+lMax*d
		if(nextPosition < 1)
			break
		elseif(nextPosition > dataSize)
			break
    end
		deltaCurrent = leading_zeros(current^mortonCodes[i+lMax*d])
	end

	l = 0
  t = lMax
  j = 0
  doit = true
	while(t>=1 || doit)
    doit = false
		nextPosition = 1
		nextPosition = i+(l+t)*d
		if (nextPosition<1)
			nextPosition = 1
		end
		if (nextPosition > dataSize)
			nextPosition = dataSize
		end
		if(leading_zeros(current$mortonCodes[nextPosition])>deltaMin)
			l = l + t
    end
		t=Int(floor(t/2))
  end
	j = i+l*d
	j = (j<1)?1:j
	j = (j>dataSize)?dataSize:j
	if(i>j)
		first = j
		last = i
	else
		first = i
		last = j
	end
	return [first,last]
end


function findSplit(mortonCodes,rng)

  firstCode = mortonCodes[rng[1]]
  lastCode = mortonCodes[rng[2]]

  if (firstCode == lastCode)
    return (rng[1] + rng[2]) >> 1
  end

  commonPrefix = leading_zeros(firstCode $ lastCode)

  split = rng[1]
  step = rng[2] - rng[1]

  doit = true

  while(step>1 || doit)
    doit = false
    step = (step + 1) >> 1
    newSplit = split + step

    if (newSplit < rng[2])
      splitCode = mortonCodes[newSplit]
      splitPrefix = leading_zeros(firstCode $ splitCode)
      if (splitPrefix > commonPrefix)
        split = newSplit
      end
    end
  end
  return split
end



function assign_inner_nodes(tree,IDs,mortonCodes)

  dataSize = length(mortonCodes)

  for i = 1:length(IDs)-1

		if(i>1)
			rng = determine_range(mortonCodes,i)
		else
			rng = [1, dataSize]
			tree.parent[i] = -1
		end

		split = findSplit(mortonCodes,rng)

    ancestors = Array(Int,2)
		if (rng[1]==split)
			ancestors[1] = IDs[split]+dataSize-1
		else
			ancestors[1] = split
		end
		if (rng[2]==split+1)
			ancestors[2] = IDs[split+1]+dataSize-1
		else
			ancestors[2] = split+1
		end

    # println("Node $(i) $(bin(mortonCodes[i])) cvg $(rng[1])-$(rng[2]) split $(split) lc $(ancestors[1]) rc $(ancestors[2])")

    tree.chld[i] = ancestors
    tree.coverage[i] = rng
    tree.parent[ancestors[1]] = i
    tree.parent[ancestors[2]] = i

    # if (mortonCodes[i]==354906307)||(mortonCodes[i]==397834323)
      println("$(i) $((mortonCodes[i])) $(rng) $(ancestors)")
    # end

  end

  return tree
end


function getMBR(tree,objects,index)
	if(index>tree.numberOfInternalNodes)
    return objects[index-tree.numberOfInternalNodes]
	else
    if isdefined(tree.mbr,index)
      return tree.mbr[index]
    else
      return -1
    end
  end
end


function getNodeMBR(tree,objects,nodes)
	vl = getMBR(tree,objects,nodes[1])
	vr = getMBR(tree,objects,nodes[2])
  if vl==-1
    return vr
  end
  if vr==-1
    return vl
  end
	v1 = minimum([vl.v1 vr.v1],2); v2 = maximum([vl.v2 vr.v2],2)
	return MBR(v1,v2)
end


function assignBV(tree,objects)
  queued = Array(Bool,tree.numberOfInternalNodes)*false

  for i = 1:length(objects)
    parent = tree.parent[tree.numberOfInternalNodes+i]

    traverse = true

    while(traverse)
			currentQ = queued[parent]
      queued[parent] = true

			if(currentQ)
        children = tree.chld[parent]

        tree.mbr[parent] = getNodeMBR(tree,objects,children)

				parent = tree.parent[parent]
				if(parent==-1)
					traverse = false
        end
			else
				traverse = false
      end
    end
  end

  return tree
end


function generate_hierarchy(tree,IDs,morton_codes,objects)
  tree = assign_inner_nodes(tree,IDs,morton_codes)

  # println("$(tree.parent)")

  # for i=1:length(objects)*2-1
  #   println("Node $(i) $(bin(morton_codes[i])) cvg $(tree.coverage[i][1])-$(tree.coverage[i][2]) lc $(tree.chld[i][1]) rc $(tree.chld[i][2]) $(tree.mbr[i])")
  # end

  tree = assignBV(tree,objects)
  #
  # for i = 1:length(IDs)
  #   println("Obj $(IDs[i]+tree.numberOfInternalNodes) $(objects[IDs[i]])")
  # end

end

function create_index(objects,lims)
  dataSize = length(objects)
  tree = radixTree(dataSize,
                    dataSize-1,
                    Array(MBR,dataSize),
                    zeros(Int,dataSize),
                    Array(Array{Int},dataSize),
                    Array(Array{Int},dataSize),
                    zeros(Int,dataSize*2-1),
                    Array(MBR,dataSize-1))

  tree.objects = deepcopy(objects)

  morton_codes = get_morton_codes(objects,lims)

  IDs = sortperm(morton_codes)
  morton_codes = morton_codes[IDs]

  generate_hierarchy(tree,IDs,morton_codes,objects)

  return tree
end



function overlap(tree::radixTree,objects::Array{MBR},ind::Int,obj::MBR)

  mbr1 = getMBR(tree,objects,ind)
  mbr2 = obj

  return prod(map(x->x<=0,(mbr1.v1-mbr2.v2).*(mbr1.v2-mbr2.v1)))
end


function probe(tree::radixTree,obj::MBR)
  pairs = Array(Int,0)
  nodestack = Array(Int,0)

  objects = tree.objects

  # println(objects)

  doit = true
  currentid = 1

  counter = 1

  # while(length(nodestack)>0 || doit)
    # doit = false
  while true

    chld = tree.chld[currentid]

    overlapL = overlap(tree,objects,chld[1],obj)
    overlapR = overlap(tree,objects,chld[2],obj)


    leftLeaf = chld[1]>tree.numberOfInternalNodes
    rightLeaf = chld[2]>tree.numberOfInternalNodes


    if(overlapL&&leftLeaf)
      push!(pairs,chld[1]-tree.numberOfInternalNodes)
    end
    if(overlapR&&rightLeaf)
      push!(pairs,chld[2]-tree.numberOfInternalNodes)
    end


    traverseL = (overlapL&&!leftLeaf)
    traverseR = (overlapR&&!rightLeaf)

    if (!traverseL && !traverseR)
      if length(nodestack)==0
        break
      end
      currentid = pop!(nodestack)
    else
      currentid = (traverseL) ? chld[1]:chld[2]
      if (traverseL && traverseR)
        push!(nodestack,chld[2])
      end
    end

    # println("$(counter) $(chld) $(overlapL) $(overlapR) $(leftLeaf) $(rightLeaf) $(traverseL) $(traverseR) $(length(nodestack))")
    #
    # counter += 1
    # # println("$(counter) $(length(nodestack))")
    # if counter > 10
    #   break
    # end

  end
  # println(pairs)
  return pairs
end

end
