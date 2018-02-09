current_path = pwd()
push!(LOAD_PATH, "$(current_path)/src")

using CoverageMapProject
using JLD
using MapVis
using runSearchPath
using LocTrack
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


# K = 10
# signalrecords = []
# realpath = []
# estimatedpath = []
# error = []
# commonerror = []
# path = "/Users/NIKMC/Project/MasterThesis/innoPositioning-Julia/test_out6_2AP/clients_jld/"
# for i=1:1
#   signal, real_path, ep, err, common_error = JLD.load("$(path)client_$(K)_$(i).jld","signalrecords","realpath","estimatedpath","error","commonerror")
#   push!(signalrecords,signal)
#   push!(realpath,real_path)
#   push!(estimatedpath,ep)
#   push!(error,err)
#   push!(commonerror,common_error)
#   println(commonerror[i])
# end

# signalrecords = []
# realpath = []
# estimatedpath = []
# error = []
# commonerror = []
# push!(signalrecords,signal)
# push!(realpath,real_path)
# push!(estimatedpath,ep)
# push!(error,err)
# push!(commonerror,common_error)

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
  name = readline()
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

MapVis.visualizeWallVis(proj)

# for ap_ind=1:length(proj.APs)
#   println("ap_ind = $(ap_ind)")
#   params = CoverageMapProject.fit_parameters(proj,ap_ind)
#   CoverageMapProject.calculate_coverage_map(proj,parameters = params)
#   # params = [147.55,-20*log10(2.4e9),0.,-0.,-2.5,-12.53,-100.]
#   # @time someFunction
# end

params = []
# params = [147.55,-20*log10(2.4e9),0.,-0.,-3.,-9.51,-41.14]    #better than another but not good
# params = [147.55,-20*log10(2.4e9),10.,-0.,-3.,-9.51,-41.14]    #not good
# params = [147.55,-20*log10(2.4e9),20.,-0.,-3.,-9.51,-41.14]    # badly
# params =  [147.55,-20*log10(2.4e9),20.,-0.,-2.5,-12.53,-12.]    # very bad
# params =  [147.55,-20*log10(2.4e9),10.,-0.,-2.5,-12.53,-12.]    #bad
# params = [147.55, -187.604, 20.0, -13.0961, -2.0, -2.0, -51.0]
# params =  [147.55,-20*log10(2.4e9),10.,-12.200322140324198,-1.1518270170732752,-13.372525638973702,-12.]    #bad
# params =  [147.55,-20*log10(2.4e9),10.,-12.200322140324198,-1.1518270170732752,-13.372525638973702,-100.]    #bad

# params = [147.55,-20*log10(2.4e9),0.,-0.,-3.,-9.51,-41.14]

if params == []
    params = CoverageMapProject.fit_parameters(proj)
end

# params = [147.55, -187.604, 10.0, -3.09613, -2.0, -2.0, -21.0] #for second AP bad result
# params = [147.55, -187.604, 20.0, -13.0961, -2.0, -2.0, -21.] #Second version for second AP bad result
# Parameter candidate [147.55, -187.604, 10.0, 2.2472, -2.29454, -3.09633, -21.0] with cost 8.572547735642075 for two Ap bad result
CoverageMapProject.calculate_coverage_map(proj,parameters = params,from_dump = true)

ssm = runSearchPath.readSSM(proj)
ssm2 = JLD.load("$(proj.path_init_data)/ssm_1.jld","ssm")
K = 20
error_of_each_path = []
for i=1:100
  signals, real_path, ep, err, common_error = JLD.load("$(proj.path_init_data)/clients_jld/client_$(K)_$(i).jld","signalrecords","realpath","estimatedpath","error","commonerror")
  estim_path = []
  espa,steps = LocTrack.estimate_path_viterbi(signals, ssm, proj.plan)
  common_error = sum(sqrt.(sum(((real_path-espa).^2),2)))/K
  push!(error_of_each_path, common_error)
  push!(estim_path, espa)
  MapVis.plot_paths(ssm2', proj, [real_path], estim_path, true, true, i)
end
println("Finish viterbi")
for i=1:length(error_of_each_path)
  println(error_of_each_path[i])
end
# prob real_path for i p_rssi= -3286.786760150085 and p_transition= -0.591517440922515
# prob real_path for i cp= -3287.378277591007 and trellis= -0.8539768411132713
# prob real_path for i  trellis2= -3287.03993765033

ssm = runSearchPath.readSSM(proj)
espa,steps = LocTrack.estimate_path_viterbi(signals, ssm, proj.plan)
coords=[]
for i=1:length(steps)
  push!(coords,LocTrack.grid2coord(fold_index(steps[i],plan),plan))
end
signals, real_path, ep, err, common_error = JLD.load("$(proj.path_init_data)/clients_jld/client_20_5.jld","signalrecords","realpath","estimatedpath","error","commonerror")
real_steps=[]
for i=1:size(real_path,1)
  push!(real_steps,LocTrack.unfold_index(LocTrack.coord2grid([real_path[i,1],real_path[i,2]],proj.plan),proj.plan))
end

println("GOOD")
#   push!(estim_path, ep)
#   MapVis.plot_paths(ssm2, proj, [real_path], estim_path, true, true, i)
#   error = sqrt.(sum(((real_path-ep).^2),2))
#   common_error = sum(sqrt.(sum(((real_path-ep).^2),2)))/10
#   println("Error of $(i) itteration = $(common_error)")
# end
  #
  #
#
#
# @time path, signal = LocTrack.path_generation(LocTrack.acelerationDist, K, ssm, proj.plan, seed ) #proj.APs[length(proj.APs)][1:2,1]
#

K = 20
seed = [-18.258438320678234 30.981275475518892]
ssm = runSearchPath.readSSM(proj)
ssm2 = JLD.load("$(proj.path_init_data)/ssm_1.jld","ssm")
error_of_each_path = []
for i=1:100
  estim_path = []
  real_path = []

  # seed = [-17.59469408768137 22.592625787210093]
   path, signal = LocTrack.path_generation(LocTrack.acelerationDist, K, ssm, proj.plan, seed ) #proj.APs[length(proj.APs)][1:2,1]
   @time ep,steps = LocTrack.estimate_path_viterbi(signal, ssm, proj.plan)
   push!(estim_path, ep)
   push!(real_path, path)
   MapVis.plot_paths(ssm2', proj, real_path, estim_path, true, true, i)
   error = sqrt.(sum(((path-ep).^2),2))
   common_error = sum(sqrt.(sum(((path-ep).^2),2)))/K
   push!(error_of_each_path, common_error)
   # println("Error of $(i) itteration = $(common_error)")
   save("$(proj.path_init_data)/clients_jld/client_$(K)_$(i).jld", "signalrecords", signal,"realpath",path,"estimatedpath", ep, "error", error, "commonerror", common_error)
   # end
 end

 println("done")

 # ssms = CoverageMapProject.load_ssms(proj)
 #
 # clients = CoverageMapProject.load_clients(proj)
 #
 # CoverageMapProject.restore_paths(clients,ssms,proj)


 for i=1:length(error_of_each_path)
   println(error_of_each_path[i])
 end
 # ssm  = []
 # ssm = runSearchPath.readSSM(proj)
 # ssm2 = JLD.load("/Volumes/FlashSD/Projects/MasterTheses/Data/indoor_test1/ssm_1.jld","ssm")
 # push!(ssm,ssm2)
 # K=10
 # seed = [-17.59469408768137 22.592625787210093]
 # @time path_real, signal = LocTrack.path_generation(LocTrack.acelerationDist, K, ssm, proj.plan, seed ) #proj.APs[length(proj.APs)][1:2,1]

# for i=1:10
#   estim_path = []
#   real_path = []
#   seed = [-18.258438320678234 30.981275475518892]
#   push!(estim_path, ep)
#   push!(real_path, path)
#   MapVis.plot_paths(ssm2', proj, real_path, estim_path, true, true, i)
#   # error = sqrt.(sum(((path-ep).^2),2))
#   # common_error = sum(sqrt.(sum(((path-ep).^2),2)))/10
#
# end
# K = 20
# signalrecords = []
# realpath = []
# estimatedpath = []
# error = []
# commonerror = []
# path = "/Users/NIKMC/Project/MasterThesis/indoor_test2/clients_jld/"
# for i=1:100
#   signal, real_path, ep, err, common_error = JLD.load("$(path)client_$(K)_$(i).jld","signalrecords","realpath","estimatedpath","error","commonerror")
#   push!(signalrecords,signal)
#   push!(realpath,real_path)
#   push!(estimatedpath,ep)
#   push!(error,err)
#   push!(commonerror,common_error)
#   println(commonerror[i])
# end


# aps = runSearchPath.readAPs()
# rssi_records = LocTrack.RssiRecord[]
# rssi_records.rssi, rssi_records.ap, rssi_records.t = readingClientTracking(aps, "$(proj.path_init_data)/clients","")
# est_path = []
# rssi_rec = []
# data_folder = "$(proj.path_init_data)/clients"
# data_path_folder = "$(proj.path_init_data)/clients_jld"
# for (fileind,file) in enumerate(readdir(data_path_folder))
#   rssi, path = JLD.load("$(proj.path_init_data)/clients_jld/client_$(fileind).jld","signalrecords","estimatedpath")
#   push!(est_path, path)
#   push!(rssi_rec, rssi)
# end
ssm2 = JLD.load("$(proj.path_init_data)/ssm_2.jld","ssm")
data_folder = "$(proj.path_init_data)/clients"
est_path = []
if length(est_path)==0
  # est_path = runSearchPath.readingClientTracking(aps, data_folder , proj)
  est_path = runSearchPath.readingBasicTracking(data_folder , proj)
end
MapVis.plot_paths(ssm2',proj,est_path, 2)


rssi, path_not = JLD.load("$(proj.path_init_data)/clients_jld/client_1.jld","signalrecords","estimatedpath")
push!(rssi_rec, rssi)

for path in est_path
  for i=1:size(path,1)-1
    println("for path between 2 point x1 = $(path[i,1]) and x2=$(path[i+1,1]) | y1 = $(path[i,2]) and y2=$(path[i+1,2])")
    ds = sqrt((path[i,1]-path[i+1,1])^2 +(path[i,2]-path[i+1,2])^2 )
    println("path = ",ds)
    dt = rssi_rec[1][i+1].t - rssi_rec[1][i].t
    println("time = ",dt)
    velocity = ds/dt
    println(" common  velocity$(i)-$(i+1) = ",velocity)
    v = (path[i,1]-path[i+1,1]) / dt
    println("velocity for x = $(v)")
    println("-------------------------------------")

  end
end
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
# traj = LocTrack.estimate_path_viterbi(signalsOfRSSI, [ssm], spaceInfo)
# pths = [ [ traj[i,1],traj[i,2] ] for i=1:size(traj,1) ]
 function synthetic_test()

for i in 10
   ssm = runSearchPath.readSSM(proj)
   path, signal = LocTrack.path_generation(LocTrack.acelerationDist, 10, ssm, proj.plan, proj.APs[length(proj.APs)][1:2,1] )
   est_path,step = LocTrack.estimate_path_viterbi(signal, ssm, proj.plan)
   MapVis.plot_paths(proj, path, est_path, true, real = true, i)
 end
 end
