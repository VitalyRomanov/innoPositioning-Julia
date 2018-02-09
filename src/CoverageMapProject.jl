module CoverageMapProject


using MapPrimitives
using MapBuilder
using MapPlan
using ImageTree
using HDF5, JLD
using DataFrames
using Optim
using MapVis
using LocTrack


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

function read_measurements(project)
    measurs = Dict()
    measur_path = "$(project.path_init_data)/loc/"
    for ap_id = 1:length(project.APs)
        if isdir("$(measur_path)/$(ap_id)")
            measurs[ap_id] = read_ap_measur(project.path_init_data,ap_id)
        end
    end
    return measurs
end

# function read_measurements(path, ap_id)
#   # import empirical RSSI measurements from disk
#   measur_avail = isdir("$(path)/loc/$(ap_id)")
#   measur = Array{Measurement}(0)
#   if measur_avail
#     for file in readdir("$(path)/loc/$(ap_id)")
#       if file[end-3:end] == ".txt" && !isnull(tryparse(Int,file[1:end-4]))
#         loc_id = parse(Int,file[1:end-4])
#         location = mean(Array(readtable("$(path)/loc/$(ap_id)/$(file)",header = false)),1)
#         rssi = Array(readtable("$(path)/rssi/$(ap_id)/$(file)",header = false)[:,2])
#         push!(measur,Measurement(loc_id,location[:],rssi))
#       end
#     end
#   end
#
#   return measur_avail,measur
# end


function read_ap_measur(path, ap_id)
  # import empirical RSSI measurements from disk
  measur = Array{Measurement}(0)
    for file in readdir("$(path)/loc/$(ap_id)")
      if file[end-3:end] == ".txt" && !isnull(tryparse(Int,file[1:end-4]))
        loc_id = parse(Int,file[1:end-4])
        location = mean(Array(readtable("$(path)/loc/$(ap_id)/$(file)",header = false)),1)
        rssi = Array(readtable("$(path)/rssi/$(ap_id)/$(file)",header = false)[:,2])
        push!(measur,Measurement(loc_id,location[:],rssi))
      end
    end

  return measur
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


function calculate_coverage_map(project::CMProject;parameters = [147.55,-20*log10(2.4e9),20.,-0.,-2.5,-12.53,-12.], from_dump = false)

  for ap_ind = 1:length(project.APs)
    AP = project.APs[ap_ind]

    ssm_path = "$(project.path_init_data)/ssm_$(ap_ind).jld"
    ssm_map_path = "$(project.path_init_data)/map_$(ap_ind)"
    dump_path = "$(project.path_init_data)/ssm_dump_$(ap_ind).txt"

    ssm = nothing

    if from_dump
        if !isfile(dump_path)
            println("No dump found!")
        else
            ssm = load_from_dump(project.plan,dump_path,parameters)
        end
    end

    if !isfile(ssm_path)
        if ssm==nothing
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
        end
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


function load_from_dump(plan,dump_path,params)
    limits = plan.limits
    size1 = limits[1,2]-limits[1,1]; size2 = limits[2,2]-limits[2,1];
    ssm = zeros(Float64,size1,size2)-3080.

    dump_ssm = readtable(dump_path, header=false)

    loc = (0.,0.); locs = []; l_paths = [];

    for row_id = 1:size(dump_ssm,1)
        c_loc = (dump_ssm[row_id,1],dump_ssm[row_id,2])
        if c_loc != loc
            loc = c_loc
            push!(locs,loc)
            push!(l_paths,0)
        end
        l_paths[end] += 1
    end

    common = sum(params[1:4])
    space_par = params[5:end]

    for l_id = 1:length(locs)
        start = sum(l_paths[1:l_id])-l_paths[l_id] + 1
        p_dump = Array(dump_ssm[start:start+l_paths[l_id]-1,3:end])
        in_build = find(x->x<2,p_dump[:,3]) # filter paths that go through buildings
        if length(in_build)==0 continue end # skip if no valid paths
        p_dump = p_dump[in_build,:]
        p_dump[:,1] = 10.0*log10.(p_dump[:,1])
        power = 10*log10(sum(10.^( (p_dump*space_par + common) / 10.))  + 1.0e-308)
        loc = locs[l_id]
        ssm[loc[1],loc[2]] = power
    end

    return ssm
end



# function recalculate_coverage_map(project::CMProject,ap_ind;parameters = [147.55,-20*log10(2.4e9),20.,-0.,-2.5,-12.53,-12.])
#   project.ssms[ap_ind] = MapBuilder.caclulate_signal_strength_matrix(project.ssms[ap_ind], project.image_trees[ap_ind],project.plan,parameters)
#
#   save_proj(project)
# end

# ============ Functions for fitting parameters ==================


# function fit_parameters(project,ap_id)
#   # X stores matrices that contain information about different signal paths
#   # to a measurement point. The matrix X[i] has 7 columns (parameters) and
#   # N rows (paths to the current location)
#   # Y stores M observations of rssi at a measurement point
#   X = Array{Array{Float64}}(0)
#   Y = Array{Array{Float64}}(0)
#
#   measur_avail,measurements = read_measurements(project.path_init_data, ap_id)
#   if ~measur_avail
#       println("No measurements available")
#       return
#   end
#
#   AP = project.APs[ap_id]
#   ap_vis = MapPlan.apVisibIndex(project.plan,AP)
#   im_tree = ImageTree.buildImageTree(project.plan,
#                                         AP,
#                                         ap_vis)
#
#
#
#   dist = [] #distance from meas. point to AP
#
#   # initial value for pathloss parameters
#   # the firts several parameters are fixed
#   # see fmin() for details
#   theta = [147.55,-20*log10(2.4e9),10.,-0.,-2.5,-12.53,-12.]
#   println(" count measur =  $(length(measurements))")
#
#   for (meas_ind,measur) in enumerate(measurements)
#     println("measur =  $(meas_ind), $(measur.rssi)")
#     x = MapBuilder.signal_paths_info(measur.location[:],
#                                 im_tree,
#                                 project.plan)
#     paths = MapBuilder.calculate_paths(measur.location[:],
#                                 im_tree,
#                                 project.plan)
#     MapVis.plot_paths(project,paths,meas_ind)
#     x = [ones(Float64,size(x,1),4) x]
#     push!(X,x)
#     push!(Y,measur.rssi)
#     dist = [dist;norm(project.APs[ap_id]-measur.location)]
#   end
#
#   theta,min_val = fmin(theta,X,Y)
#   # divide cost function by the number of observations
#   rms_error = sqrt(min_val/sum(map(i->length(Y[i]),1:length(Y))))
#   println("\nRMS $(rms_error); Optimal Parameters:\nGain: $(theta[4])\nAttenuation exponent: $(theta[5])\nReflection coefficient: $(theta[6])\nTransmission coefficient: $(theta[7])\n")
#   #calculate estimated signal strength using fit parameters
#   mean_rssi = map(i->10*log10(sum(10.^(X[i]*theta/10))),1:length(X))
#
#   MapVis.plot_fit_error(project,dist,mean_rssi,ap_id,measurements)
#   println("theta = $(theta)")
#   return theta
# end

function fit_parameters(project;n_rand_inits = 1000)#ap_id
  # X stores matrices that contain information about different signal paths
  # to a measurement point. The matrix X[i] has 7 columns (parameters) and
  # N rows (paths to the current location)
  # Y stores M observations of rssi at a measurement point
  X = Array{Array{Float64}}(0) # estimated RSSI
  Y = Array{Array{Float64}}(0) # measured RSSI

  measurements = read_measurements(project)
  if length(measurements) == 0
      println("No measurements available")
      return
  end

  dist = [] #distance from meas. point to AP

  for (ap_id,measur) in measurements
    #   iterate through different AP
    # need to build an image tree to calculate analytical signal strength
      AP = project.APs[ap_id]
      ap_vis = MapPlan.apVisibIndex(project.plan,AP)
      im_tree = ImageTree.buildImageTree(project.plan,
                                            AP,
                                            ap_vis)
    # iterate through measutement positions
      for (meas_ind,m) in enumerate(measur)
          x = MapBuilder.signal_paths_info(m.location[:],
                                      im_tree,
                                      project.plan)
          paths = MapBuilder.calculate_paths(m.location[:],
                                      im_tree,
                                      project.plan)
        #   MapVis.plot_paths(project,paths,meas_ind)
          x = [ones(Float64,size(x,1),4) x]
          push!(X,x) # add path params for measurement position meas_ind to X
          push!(Y,m.rssi) # add empirical RSSI to Y
          dist = [dist;norm(project.APs[ap_id]-m.location)]
      end
  end

  # number of observations
  n_obs = sum(map(i->length(Y[i]),1:length(Y)))

  println("Trying to find optimal parameters in $(n_rand_inits) iterations")
  # initial value for pathloss parameters
  # the firts several parameters are fixed
  # see fmin() for details
  # we change only the last three parameters as they describe the environment
  lim_bounds = [
        -4 -2;
        -30 -5;
        -50 -15
  ]
  theta = [147.55,-20*log10(2.4e9),0.,-0.,-2.5,-12.53,-12.]
  # theta = [147.55,-20*log10(2.4e9),20.,-0.,-2.5,-12.53,-12.]
  best_cost = 1e10; best_par = []
  for i=1:n_rand_inits
      sample_par = rand(3).*(lim_bounds[:,2]-lim_bounds[:,1]) + lim_bounds[:,1]
      th_temp = theta;   th_temp[5:7] = sample_par
      th_temp,min_val = fmin(theta,X,Y)
      if min_val < best_cost
          best_cost = min_val
          best_par = th_temp
          println("Parameter candidate $(best_par) with cost $(sqrt(min_val/n_obs))")
      end
  end

  theta = best_par
  min_val = best_cost
  # divide cost function by the number of observations
  rms_error = sqrt(min_val/n_obs)
  println("\nRMS $(rms_error); Optimal Parameters:\nGain: $(theta[4])\nAttenuation exponent: $(theta[5])\nReflection coefficient: $(theta[6])\nTransmission coefficient: $(theta[7])\n")
  #calculate estimated signal strength using fit parameters
  mean_rssi = map(i->10*log10(sum(10.^(X[i]*theta/10))),1:length(X))

  # MapVis.plot_fit_error(project,dist,mean_rssi,ap_id,measurements)

  return theta
end


function fmin(theta,X,Y)
  # pos-1 parameters are left from the optimization procedure

  const pos = 4
  # theta = [147.55,-20*log10(2.4e9),10.,-0.,-2.5,-12.53,-12.]
  # lower = [1,1,1,-30.,-4.,-20,-100]
  # upper = [1,1,1,30.,-2.,-2,-2]
  cost(th) = sum(map(i->sum((10*log10(sum(10.^(X[i]*[theta[1:pos-1];th]/10))) - Y[i]).^2),1:length(X)))
  # cost(th) = sum(map(i->sum((10*log10(sum(10.^(X[i]*[theta[1:pos-1];th]/10))) - Y[i]).^2)/length(Y[i]),1:length(X)))
  # res = optimize(cost, theta[pos:7], lower[pos:7], upper[pos:7], Fminbox{LBFGS}())
  res = optimize(cost,theta[pos:7],LBFGS())
  return [theta[1:pos-1];Optim.minimizer(res)],Optim.minimum(res)

end
end


# function read_measur(path)
#   measur_avail = isdir("$(path)/rssi")
#   measur = Array{Measurement}(0)
#   if measur_avail
#     for file in readdir("$(path)/rssi")
#       if file[end-3:end] == ".txt" && !isnull(tryparse(Int,file[1:end-4]))
#         loc_id = parse(Int,file[1:end-4])
#         location = []
#         rssi = Array(readtable("$(path)/rssi/$(file)",header = false)[:,2])
#         push!(measur,Measurement(loc_id,location[:],rssi))
#       end
#     end
#   end
#
#   return measur_avail,measur
# end

function load_ssms(project)
    ssms = Array{Array{Float64}}(0)
    for ap_ind = 1:length(project.APs)
        push!(ssms,JLD.load("$(project.path_init_data)/ssm_$(ap_ind).jld","ssm"))
    end
    return ssms
end

function load_clients(project)
    records = []
    for file in readdir("$(project.path_init_data)/clients")
        if file[end-2:end] != "csv" continue end
        client_path = "$(project.path_init_data)/clients/$(file)"
        push!(records,(file[1:end-4],load_client(client_path)))
    end
    return records
end

function load_client(path)
    Record = LocTrack.RssiRecord
    rssi = :rssilog_rssi
    ap_id = :rssilog_antenna
    time = :rssilog_timestamp
    readings = Array{Record}(0)

    client = readtable(path)
    n_rec = size(client,1)
    if n_rec>0
        push!(readings,Record(client[1,rssi],client[1,ap_id],client[1,time]))
    end
    last_t = client[1,time]
    if n_rec>1
        for i=2:n_rec
            if client[i,time]-last_t > 0 && client[i,rssi] < -10
                push!(readings,
                    Record(client[i,rssi],client[i,ap_id],client[i,time]))
                last_t = client[i,time]
            end
        end
    end
    return readings
end


function restore_paths(clients,ssms,project)
    paths_path = "$(project.path_init_data)/paths"
    mkpath(paths_path)
    # clients = clientz[7:8]
    for (client_fn,records) in clients
        print("\rProcessign $(client_fn)")
        path = LocTrack.estimate_path_viterbi(
                    records,
                    ssms,
                    project.plan
        )
        writetable("$(paths_path)/$(client_fn).csv",
                    convert(DataFrame,path),
                    header=false)
        MapVis.plot_paths(project,[path[:,1:2]],client_fn)
    end
    println("")

end

end
