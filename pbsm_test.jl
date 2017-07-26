current_path = "/home/vromanov/dev/innoPositioning-Julia"
current_path = "/Users/LTV/dev_projects/innoPositioning-Julia"
push!(LOAD_PATH, "$(current_path)/src")

using PBSM
using MapPrimitives
using Geometry

function load_walls3D(data_file)
  walls = Array(Wall3D,0)
  input = open(data_file)
  lims = ones(Float64,3)*[1.e50 -1.e50]
  for (line_ind,line) in enumerate(eachline(input))
    strvec = split(line,",")
    v = map(x->parse(Float64,x),strvec)

    polygon = Array(Array{Float64},0)
    coords = reshape(v[:],3,4)
    for i=1:size(coords,2)
      push!(polygon,coords[:,i])
    end
    new_wall = Wall3D(line_ind,polygon,[0,0,0,0])

    MapPrimitives.get_plane_equation!(new_wall)
    push!(walls,new_wall)

    lims[:,1] = minimum([coords lims[:,1]],2)
    lims[:,2] = maximum([coords lims[:,2]],2)
  end

  lims[:,1] = floor(lims[:,1])
  lims[:,2] = ceil(lims[:,2])

  println("$(length(walls)) walls imported")
  return walls,convert(Array{Int},lims)
end


walls,lims = load_walls3D("$(current_path)/res/coverage/tables/walls.txt")

# lims[:,1] = [-60,-30,-40]

index = PBSM.create_index(PBSM.obj2mbr(walls,wall2mbr),lims)

pairs = PBSM.probe(index,Geometry.MBR([-1,-1,-1],[5,5,5]))
