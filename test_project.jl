current_path = pwd()#"/Users/LTV/dev_projects/innoPositioning-Julia"
# current_path = "/home/vromanov/dev/innoPositioning-Julia"
# current_path = "/Users/LTV/dev_projects/innoPositioning-Julia"
cd(current_path)
push!(LOAD_PATH, "$(current_path)/src")

using CoverageMapProject
using Plots
using JLD
using LocTrack

sectorSize = 30.

last_session_available = false
options = 2
if isfile("last_session.jld")
  last_session_available = true
  options = 3
end



resp = 0
while !(resp in 1:options)
  println("\n\nChoose action")
  print("\t1 - create new project\n\t2 - load existing project\n")
  if last_session_available
    print("\t3 - load last project\n")
  end
  print("Enter option: ")
  # resp = 1
  resp = parse(Int,readline())
end


if resp==1
  print("Enter new project name: ")
  name = readline()[1:end-1]
  print("Enter path for initial data: ")
  load_path = strip(readline())
  # print("Enter saving location: ")
  # save_path = strip(readline())
  # name = "town"
  # # name = "6floor"
  # # load_path = "/Users/LTV/Documents/coverage/6floor"
  # # save_path = "/Users/LTV/Documents/coverage/6floor"
  # # load_path = "/home/vromanov/Documents/coverage/6floor"
  # # save_path = "/home/vromanov/Documents/coverage/6floor"
  # load_path = "/home/vromanov/Documents/coverage/town"
  # save_path = "/home/vromanov/Documents/coverage/town"
  proj = CoverageMapProject.create_project(load_path,
                            load_path,
                            name,
                            secSize = sectorSize)
  CoverageMapProject.save_session("$(load_path)/$(name).jld")
  CoverageMapProject.calculate_image_trees(proj)
elseif resp==2
  print("Enter existing project location :")
  proj_path = strip(readline())
  #"$(current_path)/res/coverage/init2/1/test_name.jld"
  proj = CoverageMapProject.load_project(proj_path)
  CoverageMapProject.save_session(proj_path)
elseif resp==3
  proj = CoverageMapProject.load_session()
else
  println("Unknown choice")
end

# CoverageMapProject.visualizeWallVis(proj)

# proj.image_trees_ready = false
# CoverageMapProject.calculate_image_trees(proj)

# CoverageMapProject.calculate_coverage_map(proj)
# CoverageMapProject.plot_map(proj,1)

params = CoverageMapProject.fit_parameters(proj,1)
# JLD.save("/home/ltv/Documents/coverage/town/vm.jld","vm",proj.plan.vis_matr)
# proj.ssms_ready_count = 0
# params = [147.55,-20*log10(2.4e9),0.,-0.,-2.5,-12.53,-100.]
CoverageMapProject.calculate_coverage_map(proj,parameters = params)
# CoverageMapProject.recalculate_coverage_map(proj,1,parameters = params)
# CoverageMapProject.plot_map(proj,1)
signals = []




 # measur_avail,measur = CoverageMapProject.read_measur("/Users/NIKMC/Project/MasterThesis/innoPositioning-Julia/test7")

real_path = []
proj.measur
time = [1,3,5,7,10,12,17,22,25,26,31,34,39,41]
i=1

for (meas_ind,measur) in enumerate(proj.measur)
  push!(signals,[median(measur.rssi),time[i]])
  push!(real_path, [measur.location[1],measur.location[2], measur.loc_id])
  i+=1
  # println("signals=", measur.location," | ", measur.loc_id)
end
# for (meas_ind,m) in enumerate(measur)
#   push!(signals,[median(m.rssi),m.loc_id])
#   # push!(real_path, [measur.location[1],measur.location[2], measur.loc_id])
#   # println("signals=", measur.location," | ", measur.loc_id)
# end

# println("signals=", signals)


# println("real_path=", real_path)
# proj
sort!(signals, by = x -> x[2]);
# println("signals=", signals)
sort!(real_path, by = x -> x[3]);
# println("signals=", real_path)


signalsOfRSSI = LocTrack.RssiRecord[]
for (value, id) in signals
  push!(signalsOfRSSI, LocTrack.RssiRecord(value,1,id))
end
println("signalsOfRSSI=", signalsOfRSSI)
# ssm = JLD.load("ssm1....",ssm)
ssm = JLD.load("/Users/NIKMC/Project/MasterThesis/innoPositioning-Julia/test7/ssm_1.jld","ssm")


space = proj.plan.limits[1:2,:]
grid_size = 1.
grid = convert(Array{Int},floor.((space[:,2] - space[:,1]) / grid_size))
spaceInfo = LocTrack.SpaceInfo(space,grid,grid_size)

# signal_map::Array{Array{Float64}}, +
# spaceInfo::SpaceInfo;
traj = []
traj,~ = LocTrack.estimate_path_viterbi(signalsOfRSSI, [ssm], spaceInfo)
pths = [ [ traj[i,1],traj[i,2] ] for i=1:size(traj,1) ]

# CoverageMapProject.plot_paths(proj,[pths],1)
# plot(proj.ssms[1]',seriestype=:heatmap,seriescolor=ColorGradient([colorant"white", colorant"orange", colorant"red"]),zlims=(-40,30),legend = false,grid=false,axis=false)

LocTrack.path_generation( size(signalsOfRSSI),[ssm],spaceInfo)#seed,

function plot_map(ssm,project,map_ind,filename,  paths2,paths)
    rssi_min = -100.
    rssi_max = maximum(ssm)

    println("Minimum: $(rssi_min)   Maximum: $(rssi_max)")

    plot(ssm,
        seriestype=:heatmap,
        seriescolor=ColorGradient([colorant"white", colorant"orange", colorant"red"]),
        zlims=(rssi_min,rssi_max),
        legend = false,
        grid=false,
        axis=false)

    CoverageMapProject.plot_walls!(project)

    # Plot the marker for AP location
    plot!([project.APs[map_ind][1]-project.plan.limits[1,1]],
        [project.APs[map_ind][2]-project.plan.limits[2,1]],
        markershape=:diamond,
        markercolor=:pink)

    for path in paths
      for i=1:length(path)-1
        plot!([path[i][1],path[i+1][1]]-project.plan.limits[1,1],[path[i][2],path[i+1][2]]-project.plan.limits[2,1],linecolor=:green)
      end
    end
    for path in paths2
      for i=1:length(path)-1
        plot!([path[i][1],path[i+1][1]]-project.plan.limits[1,1],[path[i][2],path[i+1][2]]-project.plan.limits[2,1],linecolor=:red)
      end
    end
    for path in paths2
      for i=1:length(path)
        x = path[i][1]-project.plan.limits[1,1]
        y = path[i][2]-project.plan.limits[2,1]
        plot!([x],[y],markershape=:diamond,markercolor=:red,subplot = 1)
        annotate!(x,y+2*7/10,text("$i",7),subplot = 1)
      end
    end
    for path in paths
      for i=1:length(path)
        x = path[i][1]-project.plan.limits[1,1]
        y = path[i][2]-project.plan.limits[2,1]
        plot!([x],[y],
            markershape=:diamond,
            markercolor=:green,
            subplot = 1)
        annotate!(x,y+2*7/10,text("$i",7),subplot = 1)
      end
    end
    savefig(filename)
    # savefig("map_.png")
end

plot_map(ssm,proj,1,"/Users/NIKMC/Project/MasterThesis/innoPositioning-Julia/test7/paths/111.png",[pths],[real_path])
# plot_map(ssm,proj,1,"/Users/NIKMC/Project/MasterThesis/innoPositioning-Julia/test8/paths/111.png",[pths])#,[real_path)

# type RssiRecord
#   rssi::Float64
#   ap::Int
#   dt::Float64
# end
#
# type SpaceInfo
#   space::Array{Int}
#   grid::Array{Int}
#   grid_size::Float64
# end
