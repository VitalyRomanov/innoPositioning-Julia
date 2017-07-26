module PBSM

using Geometry

const nod = 3 # number of dimensions
const grd_scl = 30.


function obj2mbr(walls,obj2mbr)
  mbrs = Array(MBR,length(walls))
  for (wall_ind,wall) in enumerate(walls)
    mbrs[wall_ind] = obj2mbr(wall)
  end
  return mbrs
end

type Sector
  id::Int
  objects::Array{Int}#Array{MBR}
  # geometry::Array{Array{Float64}}
  mbr::MBR
end


type Pbsm
  indexDatasetSize::Int
  objects::Array{MBR}

  grid_scale::Float64
  space_size::Array{Int}
  grid_size::Array{Int}
  sectors::Array{Sector}
  lims::Array{Int}
end


function overlap(mbr1::MBR,mbr2::MBR)
  return prod(map(x->x<=0,(mbr1.v1-mbr2.v2).*(mbr1.v2-mbr2.v1)))
end


function get_sector_ind(coord::Array{Float64},index::Pbsm)
  loc = convert(Int,ceil((coord - index.lims[1:2,1])/index.grid_scale))
  return (loc[2]-1)*index.grid_size[1]+loc[1]
end

# function get_2d_sector_ind(coord,sector_size,sectors)
#   loc = convert(Int,ceil((coord - index.lims[1:2,1])/index.grid_scale))
#   return loc
# end

function grid_to_index(grid_coord::Array{Int},index::Pbsm)
  # println(grid_coord," ",[1,index.grid_size[1],index.grid_size[2]*index.grid_size[1]])
  ind = (grid_coord)'*[1,index.grid_size[1],index.grid_size[2]*index.grid_size[1]]+1
  return ind[1]
end

function coord_to_sec_index(coord::Array{Float64},index::Pbsm)
  # test with negative coordinates
  grid_coord = convert(Array{Int},floor((coord - index.lims[:,1])/index.grid_scale))
  sector_index = (grid_coord)'*[1,index.grid_size[1],index.grid_size[2]*index.grid_size[1]]+1
  return sector_index[1]
end

function coord_to_sec_coord(coord::Array{Float64},index::Pbsm)
  # println("coord",coord)
  # println("coord",coord/index.grid_scale)
  return convert(Array{Int},floor((coord-index.lims[:,1])/index.grid_scale))
end

function prepare_sectors!(index::Pbsm)

  sector_count = 1

  for z_grid = 1:index.grid_size[3]
    for y_grid = 1:index.grid_size[2]
      for x_grid = 1:index.grid_size[1]
        geometry = Array(Array{Float64},4)
        v1 = [x_grid-1,y_grid-1,z_grid-1]*index.grid_scale+index.lims[:,1]
        v2 = [x_grid,y_grid,z_grid]*index.grid_scale+index.lims[:,1]
        index.sectors[sector_count] = Sector(sector_count,Array(MBR,0),MBR(v1,v2))
        # loc = [x_grid,y_grid,z_grid]*ones(Int,4)
        # loc = loc - [1,0,0,1;1,1,0,0]
        # for (ind,element) in enumerate(geometry)
        #   element = loc[:,ind]*index.grid_scale
        # end
        # index.sectors[sector_count] = Sector(Sector,Array(Wall,0),geometry)
        # println(sector_count," ",coord_to_sec_index(v1,index))
        sector_count+=1
      end
    end
  end
end



function sector_intersected(line::Line,sector::Sector)
  for i = 1:length(sector.geometry)-1
    if lines_crossed(line,Line(sector.geometry[i],sector.geometry[i+1]))
      return true
    end
  end
  if lines_crossed(line,Line(sector.geometry[end],sector.geometry[1]))
    return true
  end
  return false
end

function find_intersected_sectors(object::MBR,index::Pbsm)
  ind = [coord_to_sec_coord(object.v1,index) coord_to_sec_coord(object.v2,index)]

  # println("Object ",object)
  # println("sapce ",ind)

  intersected_sectors = Array(Int,0)

  ind[:,1] = maximum([ind[:,1] zeros(Int,nod,1)],2)
  # println([ind[:,1] zeros(Int,nod,1)])
  # println(maximum([ind[:,1] zeros(Int,nod,1)],2),"\n\n")

  for x=ind[1,1]:ind[1,2]
    for y=ind[2,1]:ind[2,2]
      for z=ind[3,1]:ind[3,2]
        si = grid_to_index([x,y,z],index)
        # println(si)
        if overlap(index.sectors[si].mbr,object)
          append!(intersected_sectors,si)
        end
      end
    end
  end

  # for (sec_ind,sector) in enumerate(index.sectors)
  #   if overlap(sector.mbr,object)
  #     append!(intersected_sectors,sec_ind)
  #   end
  # end

  return intersected_sectors

  # for x = ind[1,:]
  #   for y = ind[2,:]
  #     sector_index = (y-1)*index.grid_size[1]+x
  #     if sector_intersected(Line(object.v1,object.v2),index.sectors[sector_index])
  #       append!(intersected_sectors,sector_index)
  #     end
  #   end
  # end
end



function place_objects_on_grid!(index::Pbsm)
  for (obj_ind,object) in enumerate(index.objects)
    sind = find_intersected_sectors(object,index)
    for sector in index.sectors[sind]
      push!(sector.objects,obj_ind)
    end
  end
end


function create_index(objects,lims)
  dataSize = length(objects)
  space_size = lims[1:nod,2]-lims[1:nod,1]

  for i=1:length(space_size)
    if space_size[i]==0
      space_size[i] = 1
    end
  end

  grid_scale = grd_scl

  grid_size = convert(Array{Int},ceil(space_size/grid_scale))

  number_of_sector = prod(grid_size)

  sectors = Array(Sector,number_of_sector)

  index = Pbsm(dataSize,
                objects,
                grid_scale,
                space_size,
                grid_size,
                sectors,
                lims)

  prepare_sectors!(index)

  place_objects_on_grid!(index)

  return index

end


# function obj_sec_coverage(index::Pbsm,obj::MBR)



function probe(index::Pbsm,obj::MBR)

  sind = find_intersected_sectors(obj,index)
  # println(sind)
  pairs = Array(Int,0)

  for sector in index.sectors[sind]
    for obj_ind in sector.objects
      # println("Try")
      # println(obj)
      # println(index.objects[obj_ind])
      # println()

      if overlap(obj,index.objects[obj_ind])
        push!(pairs,obj_ind)
      end
    end
  end

  return pairs

end




end
