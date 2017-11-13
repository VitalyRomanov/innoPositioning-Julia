module CoverageMapProject


using MapPrimitives
using MapBuilder
using MapPlan
using ImageTree
using HDF5, JLD
using DataFrames
using Optim
using MapVis


type Measurement
  # location id assiated with the measurement point
  # (matches the filename)
  loc_id::Int
  # coordinates of the measurement point
  location::Array{Float64}
  # the array of observed RSSIs at given point
  rssi::Array{Float64}
end

type CMProject
  # path that contains the initial data and where generated files are stored
  path_init_data::String
  # project name is used for one of the files
  project_name::String

  # Objects
  APs::Array{Array{Float64}} #Array of AP locations
  plan::mapPlan

end


function read_measurements(path)
  # import empirical RSSI measurements from disk
  measur_avail = isdir("$(path)/loc")
  measur = Array{Measurement}(0)
  if measur_avail
    for file in readdir("$(path)/loc")
      if file[end-3:end] == ".txt" && !isnull(tryparse(Int,file[1:end-4]))
        loc_id = parse(Int,file[1:end-4])
        location = mean(Array(readtable("$(path)/loc/$(file)",header = false)),1)
        rssi = Array(readtable("$(path)/rssi/$(file)",header = false)[:,2])
        push!(measur,Measurement(loc_id,location[:],rssi))
      end
    end
  end

  return measur_avail,measur
end



function create_project(init_path,name;secSize = 30.)
  aps = load_aps("$(init_path)/aps.txt")
  map_plan = create_map_plan(init_path,secSize = secSize)




  project = CMProject(init_path,
                      name,
                      aps,
                      map_plan)

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
  project.plan.vis_matr = MapPlan.wallVisIndex(project.plan)
end


function save_proj(project::CMProject)
  # if !isdir(project.path_init_data)
  #   mkpath(project.path_init_data)
  # end

  JLD.save("$(project.path_init_data)/$(project.project_name).jld",
                  "init_path",project.path_init_data,
                  "proj_name",project.project_name,
                  "aps",project.APs,
                  "mapplan",project.plan)
end

function load_project(path)
  path_init_data,project_name,APs,plan = JLD.load("$(path)","init_path","proj_name","aps","mapplan")
  return CMProject(path_init_data,project_name,APs,plan)
end

function create_map_plan(init_path;secSize = 30.)
  walls,lims = load_walls3D("$(init_path)/walls.txt")

  index = MapPlan.create_index(walls,lims,grid_scl = secSize)
  index2D = MapPlan.create_2d_index(walls,lims)

  plan = MapPlan.mapPlan(walls,lims,index,index2D,Array{Bool}(0))

  if !isfile("$(init_path)/vm.jld")
    # vis_matr = MapPlan.wallVisIndexProj(plan)
    vis_matr = MapPlan.wall_vis_2D(plan)
    JLD.save("$(init_path)/vm.jld","vm",vis_matr)
  else
    vis_matr = JLD.load("$(init_path)/vm.jld","vm")
  end

  plan.vis_matr = vis_matr

  return plan
end


function load_aps(path)
  aps = Array{Array{Float64}}(0)
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
  walls = Array{Wall3D}(0)
  input = open(data_file)
  lims = ones(Float64,3)*[1.e50 -1.e50]
  for (line_ind,line) in enumerate(eachline(input))
    strvec = split(line,",")
    v = map(x->parse(Float64,x),strvec)

    polygon = Array{Array{Float64}}(0)
    coords = reshape(v[1:12],3,4)
    for i=1:size(coords,2)
      push!(polygon,coords[:,i])
    end
    new_wall = Wall3D(line_ind,polygon)#,[0,0,0,0])

    # MapPrimitives.get_plane_equation!(new_wall)
    push!(walls,new_wall)

    lims[:,1] = minimum([coords lims[:,1]],2)
    lims[:,2] = maximum([coords lims[:,2]],2)
  end

  lims[:,1] = floor.(lims[:,1])
  lims[:,2] = ceil.(lims[:,2])

  println("$(length(walls)) walls imported")
  return walls,convert(Array{Int},lims)
end


# function calculate_image_trees(project::CMProject)
#   if project.image_trees_ready
#     println("Image trees ready")
#   else
#     for (ap_ind,AP) in enumerate(project.APs)
#       push!(project.AP_visibilities,MapPlan.apVisibIndex(project.plan,AP))
#       push!(project.image_trees,ImageTree.buildImageTree(project.plan,AP,project.AP_visibilities[ap_ind]))
#       save_proj(project)
#     end
#     project.image_trees_ready = true
#     save_proj(project)
#   end
# end


function calculate_coverage_map(project::CMProject;parameters = [147.55,-20*log10(2.4e9),20.,-0.,-2.5,-12.53,-12.])

  for ap_ind = 1:length(project.APs)
    AP = project.APs[ap_ind]

    ssm_path = "$(project.path_init_data)/ssm_$(ap_ind).jld"
    ssm_map_path = "$(project.path_init_data)/map_$(ap_ind)"
    dump_path = "$(project.path_init_data)/ssm_dump_$(ap_ind).txt"

    ssm = nothing

    if !isfile(ssm_path)
      println("Creating visibility index for AP $(ap_ind)")
      ap_vis = MapPlan.apVisibIndex(project.plan,AP)

      println("Creating image tree for AP $(ap_ind)")
      im_tree = ImageTree.buildImageTree(project.plan,
                                            AP,
                                            ap_vis)
      println("Calculating coverage map for AP $(ap_ind)")
      ssm = MapBuilder.caclulate_signal_strength_matrix(im_tree,
                                                        project.plan,
                                                        parameters,
                                                        dump_path)
      JLD.save(ssm_path,"ssm",ssm)
    end

    println("Coverage map for AP $(ap_ind) is ready")

    if !isfile(ssm_map_path)
      if ssm==nothing
        ssm = JLD.load(ssm_path,"ssm")
      end
      println("Plotting coverage map for AP $(ap_ind)")
      MapVis.plot_map(ssm',project,ap_ind,ssm_map_path)
    end

    save_proj(project)
  end


end

# function recalculate_coverage_map(project::CMProject,ap_ind;parameters = [147.55,-20*log10(2.4e9),20.,-0.,-2.5,-12.53,-12.])
#   project.ssms[ap_ind] = MapBuilder.caclulate_signal_strength_matrix(project.ssms[ap_ind], project.image_trees[ap_ind],project.plan,parameters)
#
#   save_proj(project)
# end



function fit_parameters(project,ap_id)
  # X stores matrices that contain information about different signal paths
  # to a measurement point. The matrix X[i] has 7 columns (parameters) and
  # N rows (paths to the current location)
  # Y stores M observations of rssi at a measurement point
  X = Array{Array{Float64}}(0)
  Y = Array{Array{Float64}}(0)

  measur_avail,measurements = read_measurements(project.path_init_data)
  if ~measur_avail
      println("No measurements available")
      return
  end

  AP = project.APs[ap_id]
  ap_vis = MapPlan.apVisibIndex(project.plan,AP)
  im_tree = ImageTree.buildImageTree(project.plan,
                                        AP,
                                        ap_vis)



  dist = [] #distance from meas. point to AP

  # initial value for pathloss parameters
  # the firts several parameters are fixed
  # see fmin() for details
  theta = [147.55,-20*log10(2.4e9),20.,-0.,-2.5,-12.53,-12.]
  for (meas_ind,measur) in enumerate(measurements)
    x = MapBuilder.signal_paths_info(measur.location[:],
                                im_tree,
                                project.plan)
    paths = MapBuilder.calculate_paths(measur.location[:],
                                im_tree,
                                project.plan)
    MapVis.plot_paths(project,paths,meas_ind)
    x = [ones(Float64,size(x,1),4) x]
    push!(X,x)
    push!(Y,measur.rssi)
    dist = [dist;norm(project.APs[ap_id]-measur.location)]
  end

  theta,min_val = fmin(theta,X,Y)
  # divide cost function by the number of observations
  rms_error = sqrt(min_val/sum(map(i->length(Y[i]),1:length(Y))))
  println("\nRMS $(rms_error); Optimal Parameters:\nGain: $(theta[4])\nAttenuation exponent: $(-theta[5])\nReflection coefficient: $(theta[6])\nTransmission coefficient: $(theta[7])\n")
  #calculate estimated signal strength using fit parameters
  mean_rssi = map(i->10*log10(sum(10.^(X[i]*theta/10))),1:length(X))

  MapVis.plot_fit_error(project,dist,mean_rssi,ap_id,measurements)

  return theta
end




function fmin(theta,X,Y)
  # pos-1 parameters are left from the optimization procedure
  const pos = 4
  cost(th) = sum(map(i->sum((10*log10(sum(10.^(X[i]*[theta[1:pos-1];th]/10))) - Y[i]).^2),1:length(X)))
  # cost(th) = sum(map(i->sum((10*log10(sum(10.^(X[i]*[theta[1:pos-1];th]/10))) - Y[i]).^2)/length(Y[i]),1:length(X)))
  res = optimize(cost,theta[pos:7],LBFGS())
  return [theta[1:pos-1];Optim.minimizer(res)],Optim.minimum(res)
end


end
