module CoverageMapProject

# using Geometry ## remove
using MapPrimitives
using MapBuilder
using MapPlan
using ImageTree
using RadixTree
using HDF5, JLD

type CMProject
  path_init_data::String
  path_save_data::String
  project_name::String

  # Objects
  APs::Array{Array{Float64}}
  AP_visibilities::Array{Array{Bool}}
  plan::mapPlan
  image_trees::Array{Array{treeNode}}

  # Flags
  image_trees_ready::Bool

end

function create_project(init_path,save_path,name)
  aps = load_aps("$(init_path)/aps.txt")
  map_plan = create_map_plan(init_path,aps)
  image_trees = Array(Array{treeNode},0)
  AP_visibilities = Array(Array{Bool},0)
  image_trees_ready = false
  project = CMProject(init_path,
                      save_path,
                      name,
                      aps,
                      AP_visibilities,
                      map_plan,
                      image_trees,
                      image_trees_ready)

  save_proj(project)

  return project
end

function save_session(path)
  save("last_session.jld","project_path",path)
end

function load_session()
  return CoverageMapProject.load_project(load("last_session.jld","project_path"))
end


function regenerate_visibility_matrix!(project::CMProject)
  project.plan.vis_matr = MapPlan.create_wall_visibility_matrix(project.plan)
end


function save_proj(project::CMProject)
  if !isdir(project.path_save_data)
    mkpath(project.path_save_data)
  end

  JLD.save("$(project.path_save_data)/$(project.project_name).jld",
                  "init_path",project.path_init_data,
                  "save_path",project.path_save_data,
                  "proj_name",project.project_name,
                  "aps",project.APs,
                  "ap_visibilities",project.AP_visibilities,
                  "mapplan",project.plan,
                  "image_trees",project.image_trees,
                  "image_trees_ready",project.image_trees_ready)
  # save("$(project.path_save_data)/$(project.project_name).jld","CMProject",project)
end

function load_project(path)
  path_init_data,path_save_data,project_name,APs,AP_visibilities,plan,image_trees,image_trees_ready = JLD.load("$(path)","init_path","save_path","proj_name","aps","ap_visibilities","mapplan","image_trees","image_trees_ready")
  return CMProject(path_init_data,path_save_data,project_name,APs,AP_visibilities,plan,image_trees,image_trees_ready)
end

function create_map_plan(init_path,aps)
  walls,lims = load_walls3D("$(init_path)/walls.txt")

  index = RadixTree.create_index(RadixTree.obj2mbr(walls,wall2mbr),lims)

  plan = MapPlan.mapPlan(walls,lims,index,Array(Bool))

  vis_matr = MapPlan.create_wall_visibility_matrix(plan)

  plan.vis_matr = vis_matr

  return plan
end


function load_aps(path)
  aps = Array(Array{Float64},0)
  aps_file = open(path,"r")


  for (line_ind,line) in enumerate(eachline(aps_file))
    strvec = split(line,",")
    v = map(x->parse(Float64,x),strvec)
    push!(aps,v)
  end
  close(aps_file)
  return aps
end


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

  println("$(length(walls)) walls imported")
  return walls,convert(Array{Int},ceil(lims))
end


function calculate_image_trees(project::CMProject)
  if project.image_trees_ready
    println("Image trees ready")
  else
    for (ap_ind,AP) in enumerate(project.APs)
      push!(project.AP_visibilities,MapPlan.create_ap_visibility(project.plan,AP))
      push!(project.image_trees,ImageTree.build_image_tree(project.plan,AP,project.AP_visibilities[ap_ind]))
      save_proj(project)
    end
    project.image_trees_ready = true
    save_proj(project)
  end
end




end


# function load_visibility_matrix(visibility_matrix_path,plan::mapPlan)
#   if isfile(visibility_matrix_path)
#     println("Loading visibility matrix")
#     visibility_matrix = load(visibility_matrix_path,"visibility_matrix")
#   else
#     println("Visibility matrix not found. Calculating...")
#     visibility_matrix = create_wall_visibility_matrix(plan)
#     save(visibility_matrix_path,"visibility_matrix",visibility_matrix)
#   end
#   return visibility_matrix
# end
