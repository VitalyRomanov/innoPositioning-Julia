module runSearchPath

# current_path = pwd()
# push!(LOAD_PATH, "$(current_path)/src")
#
# using CoverageMapProject
# using JLD
# using MapVis
#
# if isfile("last_session.jld")
#   last_session_available = true
#   options = 3
# end
# resp = 0
# while !(resp in 1:options)
#   if last_session_available
#     print("\t3 - load last project\n")
#   end
import CoverageMapProject
using DataFrames
# Measurement = CoverageMapProject.Measurement

type Measurement
  # location id assiated with the measurement point
  # (matches the filename)
  loc_id::Int
  # coordinates of the measurement point
  location::Array{Float64}
  # the array of observed RSSIs at given point
  rssi::Array{Float64}
  dt::Array{Int64}
end

function init(project)
  measurements = Array{Array{Measurement}}(0)
  real_measurements = Array{Array{Measurement}}(0)
  for iter in length(project.APs)
    measur_avail,meas = read_measurs(project.path_init_data, iter)
    real, real_measure = read_real_path(project.path_init_data,iter)
      # measur_avail,meas = read_measur(project.path_init_data)
    if ~measur_avail
        println("No measurements available $(iter)")
        return
      end
      if ~real
          println("No real measurements available $(iter)")
          return
        end
    push!(real_measurements,real_measure)
    push!(measurements,meas)
  end

for ap_ind= 1:length(project.APs)
  AP = project.APs[ap_ind]
  est_path = calculate(project,ap_ind,measurements[ap_ind])
  MapVis.plot_paths(project,real_measurements[ap_ind],est_path,true)
end

  #read_real_path(project)
end

  function calculate(project, ap_ind, measures)
    ssm_path = "$(project.path_init_data)/ssm_$(ap_ind).jld"
    space = project.plan.limits[1:2,:]
    grid_size = 1.
    #???? measur.location[] - 117
    for (meas_ind,measur) in enumerate(measures)
      println(measur.location)
      println("$(measur.rssi), $(measur.dt)")
      dt = (meas_ind==1)?1:measures[meas_ind].dt-measures[meas_ind-1].dt # RssiRecord stores dt instead of t
      iloc = Int.(ceil.((measur.location[1:2]-space[:,1])./grid_size)) # the index of ssm that corresponds to the true location
      ss = ssm[iloc[1],iloc[2]] # use signal strength from previous measurements (for now)
      # push!(signals,[measur.rssi[1],dt]) # these values are way way off
      push!(signals,[mode(measur.rssi),dt,measur.loc_id]) # use the same exact values as we expect from ssm for testing
      # push!(real_path, [measur.location[1],measur.location[2], measur.loc_id])
      # println("Expected $(ss), observed $(signals[end][1]) $(measur.loc_id)")
    end
    sort!(signals, by = x -> x[3]);
    signalsOfRSSI = LocTrack.RssiRecord[]
    for (value, id, loc_id) in signals
      println("rssi = $(value) \t loc = $(loc_id)")
      push!(signalsOfRSSI, LocTrack.RssiRecord(value,ap_ind,id))
    end
    traj = []
    traj,~ = LocTrack.estimate_path_viterbi(signalsOfRSSI, [ssm_path], getSpaceInfo(project.plan.limits[1:2,:]))
    return pths = [ [ traj[i,1],traj[i,2] ] for i=1:size(traj,1) ]
  end



  function read_real_path(path, ap_id)
    # app_id = length(project.APs)
      measur_avail = isdir("$(path)/real_path/$(ap_id)")
      measur = Array{Measurement}(0)
      if measur_avail
        for file in readdir("$(path)/real_path/$(ap_id)")
          if file[end-3:end] == ".txt" && !isnull(tryparse(Int,file[1:end-4]))
            loc_id = parse(Int,file[1:end-4])
            location =  mean(Array(readtable("$(path)/loc/$(ap_id)/$(file)",header = false)),1)
            rssi = [] #Array(readtable("$(path)/real_path/$(iter)/$(file)",header = false)[:,2])
            dt = 1
            push!(measur,Measurement(loc_id,location[:],rssi,dt))
          end
        end
      end
    return measur_avail,measur
  end

  function read_measurs(path, ap_id)
    # import empirical RSSI measurements from disk
      measur_avail = isdir("$(path)/rssi_path/$(ap_id)")
      measur = Array{Measurement}(0)            #Как тут сделать массмв массивов
      if measur_avail
        for file in readdir("$(path)/rssi_path/$(ap_id)")
          if file[end-3:end] == ".txt" && !isnull(tryparse(Int,file[1:end-4]))
            loc_id = parse(Int,file[1:end-4])
            location = []
            rssi = Array(readtable("$(path)/rssi_path/$(ap_id)/$(file)",header = false)[:,2])
            dt = Array(readtable("$(path)/rssi_path/$(ap_id)/$(file)",header = false)[:,3])
            push!(measur,Measurement(loc_id,location[:],rssi,dt))      # как сделать добавление в массив массива
          end
        end
      end
    return measur_avail,measur
  end
end
