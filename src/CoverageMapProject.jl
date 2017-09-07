module CoverageMapProject

# using Geometry ## remove
using MapPrimitives
using MapBuilder
using MapPlan
using ImageTree
using Plots
using StatPlots
# using RadixTree
using HDF5, JLD
using DataFrames
using Optim

gr()


type Measurement
  loc_id::Int
  location::Array{Float64}
  rssi::Array{Float64}
end

type CMProject
  path_init_data::String
  path_save_data::String
  project_name::String

  # Objects
  APs::Array{Array{Float64}}
  AP_visibilities::Array{Array{Bool}}
  plan::mapPlan
  image_trees::Array{Array{treeNode}}
  ssms::Array{Array{Float64}}

  # Flags
  image_trees_ready::Bool
  coverage_maps_ready::Bool
  measur_avail::Bool

  # Counters
  ssms_ready_count::Int # remove

  # Measurements
  measur::Array{Measurement}

end


function read_measurements(path)
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



function create_project(init_path,save_path,name)
  aps = load_aps("$(init_path)/aps.txt")
  map_plan = create_map_plan(init_path)
  image_trees = Array{Array{treeNode}}(0)
  AP_visibilities = Array{Array{Bool}}(0)
  ssms = Array{Array{Float64}}(0)
  image_trees_ready = false
  coverage_maps_ready = false
  ssms_ready_count = 0

  # for i=1:length(aps)
  #   push!(ssms,zeros(Float64,map_plan.limits[1,2]-map_plan.limits[1,1],map_plan.limits[2,2]-map_plan.limits[2,1]))
  # end

  measur_avail,measur = read_measurements("$(init_path)")


  project = CMProject(init_path,
                      save_path,
                      name,
                      aps,
                      AP_visibilities,
                      map_plan,
                      image_trees,
                      ssms,
                      image_trees_ready,
                      coverage_maps_ready,
                      measur_avail,
                      ssms_ready_count,
                      measur)

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
                  "ssms",project.ssms,
                  "image_trees_ready",project.image_trees_ready,
                  "coverage_maps_ready",project.coverage_maps_ready,
                  "measur_avail",project.measur_avail,
                  "ssms_ready_count",project.ssms_ready_count,
                  "measur",project.measur)
  # save("$(project.path_save_data)/$(project.project_name).jld","CMProject",project)
end

function load_project(path)
  path_init_data,path_save_data,project_name,APs,AP_visibilities,plan,image_trees,ssms,image_trees_ready,coverage_maps_ready,measur_avail,ssms_ready_count,measur = JLD.load("$(path)","init_path","save_path","proj_name","aps","ap_visibilities","mapplan","image_trees","ssms","image_trees_ready","coverage_maps_ready","measur_avail","ssms_ready_count","measur")
  return CMProject(path_init_data,path_save_data,project_name,APs,AP_visibilities,plan,image_trees,ssms,image_trees_ready,coverage_maps_ready,measur_avail,ssms_ready_count,measur)
end

function create_map_plan(init_path)
  walls,lims = load_walls3D("$(init_path)/walls.txt")

  # index = RadixTree.create_index(RadixTree.obj2mbr(walls,wall2mbr),lims)
  index = MapPlan.create_index(walls,lims,30.)

  plan = MapPlan.mapPlan(walls,lims,index,Array{Bool}(0))

  if !isfile("$(init_path)/vm.jld")
    vis_matr = MapPlan.create_wall_visibility_matrix(plan)
    JLD.save("$(init_path)/vm.jld","vm",vis_matr)
  else
    vis_matr = JLD.load("$(init_path)/vm.jld","vm")
  end

  # Profile.print(combine = true,sortedby=:count)

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
    coords = reshape(v[:],3,4)
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



# function calculate_image_tree(project::CMProject,AP)
#     ap_vis = MapPlan.create_ap_visibility(project.plan,AP);
#     im_tree = ImageTree.build_image_tree(project.plan,
#                                         AP,
#                                         ap_vis)




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


function calculate_coverage_map(project::CMProject;parameters = [147.55,-20*log10(2.4e9),20.,-0.,-2.5,-12.53,-12.])
  if project.coverage_maps_ready
    println("Coverage maps are calculated")
    return
  end

  for ap_ind = project.ssms_ready_count+1:length(project.APs)
    AP = project.APs[ap_ind]

    ssm_path = "$(project.path_save_data)/ssm_$(ap_ind).jld"
    ssm_map_path = "$(project.path_save_data)/map_$(ap_ind)"

    ssm = nothing

    if !isfile(ssm_path)
      println("Creating visibility index for AP $(ap_ind)")
      ap_vis = MapPlan.create_ap_visibility(project.plan,AP)

      println("Creating image tree for AP $(ap_ind)")
      im_tree = ImageTree.build_image_tree(project.plan,
                                            AP,
                                            ap_vis)
      println("Calculating coverage map for AP $(ap_ind)")
      ssm = MapBuilder.caclulate_signal_strength_matrix(im_tree,
                                                        project.plan,
                                                        parameters)
      JLD.save(ssm_path,"ssm",ssm)
    end

    println("Coverage map for AP $(ap_ind) is ready")

    if !isfile(ssm_map_path)
      if ssm==nothing
        ssm = JLD.load(ssm_path,"ssm")
      end
      println("Plotting coverage map for AP $(ap_ind)")
      plot_map(ssm',project,ap_ind,ssm_map_path)
    end
    # project.ssms[ap_ind] = MapBuilder.caclulate_signal_strength_matrix(project.ssms[ap_ind],
    #                                                                   project.image_trees[ap_ind],
    #                                                                   project.plan,
    #                                                                   parameters)


    # push!(project.ssms,ssm)
    project.ssms_ready_count += 1
    save_proj(project)
  end

  project.coverage_maps_ready = true
end

function recalculate_coverage_map(project::CMProject,ap_ind;parameters = [147.55,-20*log10(2.4e9),20.,-0.,-2.5,-12.53,-12.])
  project.ssms[ap_ind] = MapBuilder.caclulate_signal_strength_matrix(project.ssms[ap_ind], project.image_trees[ap_ind],project.plan,parameters)

  save_proj(project)
end


function plot_walls!(project::CMProject)
  for wall in project.plan.walls
    if wall.polygon[1][3]==wall.polygon[3][3]
      continue
    end
    xs = [wall.polygon[1][1],wall.polygon[4][1]]
    ys = [wall.polygon[1][2],wall.polygon[4][2]]
    plot!(xs,ys,
        linecolor=:black,
        xlims = project.plan.limits[1,:],
        ylims = project.plan.limits[2,:],
        legend = false)
  end
end


function plot_map(ssm,project,map_ind,filename)
    rssi_min = -100.
    rssi_max = maximum(ssm)


    # for x = 1:size(ssm,1), y = 1:size(ssm,2)
    #   if ssm[x,y] < rssi_min && ssm[x,y] > -900.
    #     rssi_min = ssm[x,y]
    #   end
    # end

    println("Minimum: $(rssi_min)   Maximum: $(rssi_max)")

    plot(ssm,
        seriestype=:heatmap,
        seriescolor=ColorGradient([colorant"white", colorant"orange", colorant"red"]),
        zlims=(rssi_min,rssi_max),
        legend = false,
        grid=false,
        axis=false)

    plot_walls!(project)

    plot!([project.APs[map_ind][1]],[project.APs[map_ind][2]],markershape=:diamond,markercolor=:pink)

    savefig(filename)
    # savefig("map_.png")
end

function plot_map(project::CMProject,map_ind)
  ssm = project.ssms[map_ind]'
  ssm_map_path = "$(project.path_save_data)/map_$(map_ind).png"

  plot_map(ssm,project,map_ind,filename)

  # rssi_min = -100.
  # rssi_max = maximum(ssm)
  #
  # # for x = 1:size(ssm,1), y = 1:size(ssm,2)
  # #   if ssm[x,y] < rssi_min && ssm[x,y] > -900.
  # #     rssi_min = ssm[x,y]
  # #   end
  # # end
  #
  # println("Minimum: $(rssi_min)   Maximum: $(rssi_max)")
  #
  # plot(ssm,
  #     seriestype=:heatmap,
  #     seriescolor=ColorGradient([colorant"white", colorant"orange", colorant"red"]),
  #     zlims=(rssi_min,rssi_max),
  #     legend = false,
  #     grid=false,
  #     axis=false)
  #
  # plot_walls!(project)
  #
  # plot!([project.APs[map_ind][1]],[project.APs[map_ind][2]],markershape=:diamond,markercolor=:pink)
  #
  # savefig("$(project.path_save_data)/map_$(map_ind).png")
end


function plot_paths(project::CMProject,paths,mp_ind)
  plot(size(600,600),
      axis = false)
  plot_walls!(project)

  for path in paths
    for i=1:length(path)-1
      plot!([path[i][1],path[i+1][1]],[path[i][2],path[i+1][2]],linecolor=:red)
    end
  end

  if !isdir("$(project.path_save_data)/paths")
    mkdir("$(project.path_save_data)/paths")
  end

  savefig("$(project.path_save_data)/paths/$(mp_ind).svg")
end


function fit_parameters(project::CMProject,ap_id::Int)
  # X stores matrices that contain information about different signal paths
  # to a measurement point. The matrix X[i] has 7 columns (parameters) and
  # N rows (paths to the current location)
  # Y stores M observations of rssi at a measurement point
  X = Array{Array{Float64}}(0)
  Y = Array{Array{Float64}}(0)

  fs = 7 # font size for plot annotations

  dist = [] #distance from meas. point to AP

  # initial value for pathloss parameters
  # the firts several parameters are fixed
  # see fmin() for details
  theta = [147.55,-20*log10(2.4e9),20.,-0.,-2.5,-12.53,-12.]
  for (meas_ind,measur) in enumerate(project.measur)
    x = MapBuilder.signal_paths_info(measur.location[:],
                                project.image_trees[ap_id],
                                project.plan)
    paths = MapBuilder.calculate_paths(measur.location[:],
                                project.image_trees[ap_id],
                                project.plan)
    plot_paths(project,paths,meas_ind)
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

  # create a plot with specified parameters
  # this plot will contain the area map in the first subplot and
  # comparison of the observed and estimated rssi values in the second
  plot(legend = false,
      grid=false,
      # axis=false,
      size=(900,1800),
      layout=grid(2,1,heights=[.5,.5]))
  plot_walls!(project)


  order = sortperm(dist)

  plot!(dist[order],mean_rssi[order],legend=false,subplot=2)
  plot!([project.APs[ap_id][1]],[project.APs[ap_id][2]],
      markershape=:diamond,
      markercolor=:pink,
      subplot = 1)

  for elem in order
    x = project.measur[elem].location[1]
    y = project.measur[elem].location[2]
    plot!([x],[y],
        markershape=:diamond,
        markercolor=:blue,
        subplot = 1)
    annotate!(x,y+2*fs/10,text("$elem",fs),subplot = 1)
    boxplot!([dist[elem]],project.measur[elem].rssi,lab="Loc $(elem)",subplot = 2)
    annotate!(dist[elem]*1.01,mean(project.measur[elem].rssi)*1.01,text("$elem"),subplot = 2)
  end

  savefig("$(project.path_save_data)/all_it_takes.svg")

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
