module PBSM

using Geometry

const nod = 2 # number of dimensions

export MBR

function obj2mbr(walls,obj2mbr)
  mbrs = Array(MBR,length(walls))
  for (wall_ind,wall) in enumerate(walls)
    mbrs[wall_ind] = obj2mbr(wall)
  end
  return mbrs
end

type Sector
  id::Int
  objects::Array{MBR}
  geometry::Array{Array{Float64}}
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


function get_sector_ind(coord::Array{Float64},index::Pbsm)
  loc = convert(Int,ceil((coord - index.lims[1:2,1])/index.grid_scale))
  return (loc[2]-1)*index.grid_size[1]+loc[1]
end

function get_2d_sector_ind(coord,sector_size,sectors)
  loc = convert(Int,ceil((coord - index.lims[1:2,1])/index.grid_scale))
  return loc
end


function prepare_sectors!(index::Pbsm)

  sector_count = 1

  for y_grid = 1:grid_size[2]
    for x_grid = 1:grid_size[1]
      geometry = Array(Array{Float64},4)
      loc = [x_grid,y_grid]*ones(Int,4)
      loc = loc - [1,0,0,1;1,1,0,0]
      for (ind,element) in enumerate(geometry)
        element = loc[:,ind]*index.grid_scale
      end
      index.sectors[sector_count] = Sector(Sector,Array(Wall,0),geometry)

      sector_count+=1
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
  ind = [get_2d_sector_ind(object,v1) get_2d_sector_ind(object,v2)]

  intersected_sectors = Array(Int,0)

  for x = ind[1,:]
    for y = ind[2,:]
      sector_index = (y-1)*index.grid_size[1]+x
      if sector_intersected(Line(object.v1,object.v2),index.sectors[sector_index])
        append!(intersected_sectors,sector_index)
      end
    end
  end



function place_objects_on_grid(index::Pbsm)
  for object in index.objects
    ind = find_intersected_sectors(object,index)
    for sector in index.sectors[ind]
      append!(sector.objects,object)
    end
  end
end


function create_index(objects,lims)
  dataSize = length(objects)
  space_size = lims[1:2,2]-lims[1:2,1]

  grid_scale = 30.

  grid_size = convert(Array{Int},ceil(space_size/grid_scale))

  number_of_sector = prod(grid_scale)

  index = Pbsm(dataSize,
                objects,
                grid_scale,
                space_size,
                grid_size,
                Array(Sector,number_of_sector),
                lims[1:2,1:2])

  prepare_sectors!(index)

  place_objects_on_grid!(index)

end


function probe(index::Pbsm,obj::MBR)

end




end
