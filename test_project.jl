current_path = pwd()
push!(LOAD_PATH, "$(current_path)/src")

using CoverageMapProject
using JLD
using MapVis
using runSearchPath
# defines the grid size for PBSM index
sectorSize = 30.

# check if have a recently launched project
# add an option to load the last project if true
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
  try
      resp = parse(Int,readline())
  catch
      println("Enter an integer number")
  end
end


if resp==1
  print("Enter new project name: ")
  name = readline()[1:end]
  print("Enter path for the project: ")
  load_path = strip(readline())
  # name = "test8"
  # load_path = "/Users/LTV/Dropbox (Innopolis)/Work/innoPositioning-Julia-dev/test8"
  proj = CoverageMapProject.create_project(load_path,
                            name,
                            secSize = sectorSize)
  CoverageMapProject.save_session("$(load_path)/$(name).jld")
elseif resp==2
  print("Enter existing project location :")
  proj_path = strip(readline())
  proj = CoverageMapProject.load_project(proj_path)
  CoverageMapProject.save_session(proj_path)
elseif resp==3
  proj = CoverageMapProject.load_session()
# elseif resp==4
  # print("Enter existing first project location :")
  # proj_path = strip(readline())
  # proj = CoverageMapProject.load_project(proj_path)
  # CoverageMapProject.save_session(proj_path)
  # print("Enter existing second project location :")
  # proj_path_2 = strip(readline())
  # proj_2 = CoverageMapProject.load_project(proj_path_2)
  # CoverageMapProject.save_session(proj_path_2)
  # params_1 = CoverageMapProject.fit_parameters(proj,1)
  # params_2 = CoverageMapProject.fit_parameters(proj_2,2)
  # CoverageMapProject.calculate_coverage_map(proj,parameters = params)
else
  println("Unknown choice")
end

# MapVis.visualizeWallVis(proj)

# for ap_ind=1:length(proj.APs)
#   println("ap_ind = $(ap_ind)")
#   params = CoverageMapProject.fit_parameters(proj,ap_ind)
#   CoverageMapProject.calculate_coverage_map(proj,parameters = params)
#   # params = [147.55,-20*log10(2.4e9),0.,-0.,-2.5,-12.53,-100.]
#   # @time someFunction
# end

# params = []
# params = [147.55,-20*log10(2.4e9),0.,-0.,-3.,-9.51,-41.14]
# if params == []
#     params = CoverageMapProject.fit_parameters(proj,1)
# end
# CoverageMapProject.calculate_coverage_map(proj,parameters = params,from_dump = true)


aps = runSearchPath.readAPs()
# rssi_records = LocTrack.RssiRecord[]
# rssi_records.rssi, rssi_records.ap, rssi_records.t = readingClientTracking(aps, "$(proj.path_init_data)/clients","")
est_path = []
rssi_rec = []
data_folder = "$(proj.path_init_data)/clients"
data_path_folder = "$(proj.path_init_data)/clients_jld"
for (fileind,file) in enumerate(readdir(data_path_folder))
  rssi, path = JLD.load("$(proj.path_init_data)/clients_jld/client_$(fileind).jld","signalrecords","estimatedpath")
  push!(est_path, path)
  push!(rssi_rec, rssi)
end


if length(est_path)==0
  est_path = runSearchPath.readingClientTracking(aps, data_folder , proj)
end
MapVis.plot_paths(proj,[est_path])
# MapVis.plot_paths(proj,real_measurements[2],[est_path],true,false)
println("Good")
# space = project.plan.limits[1:2,:]
# grid_size = 1.
# grid = convert(Array{Int},floor.((space[:,2] - space[:,1]) / grid_size))
# SpaceInfo(space,grid,grid_size)


# signals = []
#
# sm = JLD.load("/Users/NIKMC/Project/MasterThesis/innoPositioning-Julia/test10/ssm_1.jld","ssm") # use data from test 6
#
# real_path = []
# time = [1,3,5,7,10,12,17,22,25,26,31,34,39,41]
# i=1
#
# for (meas_ind,measur) in enumerate(proj.measur)
#   dt = (i==1)?1:time[i]-time[i-1] # RssiRecord stores dt instead of t
#   iloc = Int.(ceil.((measur.location[1:2]-space[:,1])./grid_size)) # the index of ssm that corresponds to the true location
#   ss = ssm[iloc[1],iloc[2]] # use signal strength from previous measurements (for now)
#   # push!(signals,[measur.rssi[1],dt]) # these values are way way off
#   push!(signals,[mode(measur.rssi),dt,measur.loc_id]) # use the same exact values as we expect from ssm for testing
#   push!(real_path, [measur.location[1],measur.location[2], measur.loc_id])
#   println("Expected $(ss), observed $(signals[end][1]) $(measur.loc_id)")
#   i+=1
# end
#
# # proj
# sort!(signals, by = x -> x[3]);
# # println("signals=", signals)
# sort!(real_path, by = x -> x[3]);
# # println("signals=", real_path)
# for (m, m2 ) in enumerate(signals)
#   println("observed $(signals[m])")
# end
#
# signalsOfRSSI = LocTrack.RssiRecord[]
# for (value, id, loc_id) in signals
#   println("rssi = $(value) \t loc = $(loc_id)")
#   push!(signalsOfRSSI, LocTrack.RssiRecord(value,1,id))
# end
# println("signalsOfRSSI=", signalsOfRSSI)
#
#
#
# traj = []
# traj,~ = LocTrack.estimate_path_viterbi(signalsOfRSSI, [ssm], spaceInfo)
# pths = [ [ traj[i,1],traj[i,2] ] for i=1:size(traj,1) ]
