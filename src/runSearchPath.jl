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
using JLD
using LocTrack
# Measurement = CoverageMapProject.Measurement

type Measurement
  loc_id::Int
  rssi::Array{Float64}
  t::Array{Int64}
end

function init(project)
  measurements = Array{Array{Measurement}}(0)
  real_measurements = Array{Array{Measurement}}(0)
  for iter=1:length(project.APs)
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

# for ap_ind= 1:length(project.APs)
  # AP = project.APs[ap_ind]
  # println(measurements)
  est_path = calculate(project,measurements)
  # MapVis.plot_paths(project,real_measurements[ap_ind],est_path,true)
# end
  return est_path, measurements, real_measurements
  #read_real_path(project)
end

  function calculate(project, measures)
    ssm_path::Array{Array{Float64}} = []
    for ap_ind=1:length(project.APs)
        ssm_path_i = "$(project.path_init_data)/ssm_$(ap_ind).jld"
        ssm = JLD.load(ssm_path_i,"ssm")
        push!(ssm_path,ssm)
    end
    signals = []
    for ap_ind=1:length(project.APs)
      for (meas_ind,measur) in enumerate(measures[ap_ind])
        println("$(measur.rssi), $(measur.t)")
        push!(signals,[mode(measur.rssi),convert(Float64,measur.t[1]), measur.loc_id, ap_ind]) # use the same exact values as we expect from ssm for testing
      end
    end
    sort!(signals, by = x -> x[2]);
    signalsOfRSSI = LocTrack.RssiRecord[]
    for (value, t, loc_id, ap_id) in signals
      println("rssi = $(value) \t loc = $(loc_id) \t ap= $(ap_id) \t  time = $(t)")
      push!(signalsOfRSSI, LocTrack.RssiRecord(value,ap_id,convert(Float64, t[1])))
    end
    traj = []
    traj = estimate_path_viterbi(signalsOfRSSI, ssm_path, project.plan )
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
            # location = []
            rssi = Array(readtable("$(path)/rssi_path/$(ap_id)/$(file)",header = false)[:,2])
            t = Array(readtable("$(path)/rssi_path/$(ap_id)/$(file)",header = false)[:,3])
            push!(measur,Measurement(loc_id,rssi,t))      # как сделать добавление в массив массива
          end
        end
      end
    return measur_avail,measur
  end

  function readAPs(path_file_ap)
    aps = Dict()
    i=1
    # for l in "/Users/NIKMC/Project/MasterThesis/ap_i.txt"|>eachline
    for l in path_file_ap|>eachline
      fields = split(l, ',')
      for ii=1:length(fields)
        aps[fields[ii]] = i
      end
    i+=1
    end
  return aps
  end

  function readSSM(project)
    ssm_path::Array{Array{Float64}} = []
    for ap_ind=1:length(project.APs)
        ssm_path_i = "$(project.path_init_data)/ssm_$(ap_ind).jld"
        ssm = JLD.load(ssm_path_i,"ssm")
        push!(ssm_path,ssm)
    end
    return ssm_path
  end


  function readingClientTracking(data_folder,project,path_file_ap)
    aps = readAPs(path_file_ap)
    ssm = readSSM(project)
    estim_path = []
    for (fileind,file) in enumerate(readdir(data_folder))
      rssi_records = LocTrack.RssiRecord[]
      println("Reading file $(file)")
      client = readtable("$(data_folder)/$(file)")
      println("Formatting...")
      rssi, ip, timestamps = get_records(client)
      ap = convertIpInID(ip,aps)
      println("Converting...")
      t = convertTimeStamp(timestamps)
      println("Checking data...")

      if(length(rssi)==length(t) && length(rssi)==length(ap) && length(t)==length(ap))
        println("it's ok...")
        for i=1:length(rssi)
          push!(rssi_records, LocTrack.RssiRecord(rssi[i],ap[i],t[i]))
        end
      end
      # push!(rssi_records, LocTrack.RssiRecord(rssi,ap,t))
      # return rssi_records, time, ap
          # rssi_records
    println("rssi_records = ", rssi_records)
    println("start estimate path with viterbi")
    @time ep = estimate_path_viterbi(rssi_records, ssm, project.plan )
    # println("end estimate path with viterbi")
    # println("path = ", ep)
    push!(estim_path, ep)

    println("Saving for client $(file)")
    save("$(project.path_init_data)/clients_jld/client_$(file).jld", "signalrecords", rssi_records,"estimatedpath",ep)
    # save("$(current_path)/res/paths/client_$(fileind).jld", "signalrecords", signal_records,"estimatedpath",ep,"timestamps", time)
  end
  # Array(readtable("$(path)/rssi_path/$(ap_id)/$(file)",header = false)[:,2])
  # timestamps = Array(client[:,5])
  # timestamps = timestamps[1:11922]
  # save("$(current_path)/res/paths/client_1.jld", "signalrecords", signal_records,"estimatedpath",ep,"timestamps", timestamps)


  # time(Libc.strptime("%Y-%m-%d %T","2015-09-18 01:01:23.24963"))
return estim_path, ssm
end

function loadingSyntheticTestPathFor1AP(proj,K)
  ssm = readSSM(proj)
  ssm2 = JLD.load("$(proj.path_init_data)/ssm_1.jld","ssm")
  error_of_each_path = []
  for i=1:100
    signals, real_path, ep, err, common_error = JLD.load("$(proj.path_init_data)/clients_jld/client_$(K)_$(i).jld","signalrecords","realpath","estimatedpath","error","commonerror")
    estim_path = []
    println("start viterbi $(i)")
    @time espa = estimate_path_viterbi(signals, ssm, proj.plan)#, seed = real_path[1,:])
    # for j=1:K
    #   println("signals in real path = $(signals[j].rssi)")
    #   println("signals in estimated path = $(ssm[1][LocTrack.coord2grid(espa[j,1:2],proj.plan)[1],LocTrack.coord2grid(espa[j,1:2],proj.plan)[2]])")
    # end
    common_error = sum(sqrt.(sum(((real_path-espa[:,1:2]).^2),2)))/K
    push!(error_of_each_path, common_error)
    push!(estim_path, espa)
    MapVis.plot_paths(ssm2', proj, [real_path], estim_path, true, true, i)
  end
  for i=1:length(error_of_each_path)
    println(error_of_each_path[i])
  end
  println("Finish viterbi")
end

function creatingSyntheticTestPathFor1AP(proj,K; seed = [-18.258438320678234 30.981275475518892])
  ssm = runSearchPath.readSSM(proj)
  ssm2 = JLD.load("$(proj.path_init_data)/ssm_1.jld","ssm")
  error_of_each_path = []
  for i=1:100
    estim_path = []
    real_path = []
    # seed = [-17.59469408768137 22.592625787210093]
     path, signal = LocTrack.path_generation(LocTrack.acelerationDist, K, ssm, proj.plan, seed ) #proj.APs[length(proj.APs)][1:2,1]
     println("start viterbi $(i)")
     @time ep = estimate_path_viterbi(signal, ssm, proj.plan)
     println("Finish viterbi")
     push!(estim_path, ep)
     push!(real_path, path)
     MapVis.plot_paths(ssm2', proj, real_path, estim_path, true, true, i)
     # error = sqrt.(sum(((path-ep[:,1:2]).^2),2))
     common_error = sum(sqrt.(sum(((path-ep[:,1:2]).^2),2)))/K
     push!(error_of_each_path, common_error)
     println("Error of $(i) itteration = $(common_error)")
     println("SavingData in clients_jld/client_$(K)_$(i).jld")
     save("$(proj.path_init_data)/clients_jld/client_$(K)_$(i).jld", "signalrecords", signal,"realpath",path,"estimatedpath", ep, "error", error, "commonerror", common_error)
     # end
  end
end

function readingBasicTracking(data_folder,project)
  ssm = readSSM(project)
  estim_path = []
  for (fileind,file) in enumerate(readdir(data_folder))
    rssi_records = LocTrack.RssiRecord[]
    println("Reading file $(file)")
    client = readtable("$(data_folder)/$(file)",header = false)
    println("Formatting...")
    rssi, ap, t = getBasicRecords(client)
    # println("Converting...")
    # t = convertTimeStamp(timestamps)
    println("Checking data...")

    if (length(rssi)==length(t) && length(rssi)==length(ap) && length(t)==length(ap))
      println("it's ok...")
      for i=1:length(rssi)
        push!(rssi_records, LocTrack.RssiRecord(rssi[i],ap[i],t[i]))
      end
    end
    # push!(rssi_records, LocTrack.RssiRecord(rssi,ap,t))
    # return rssi_records, time, ap
        # rssi_records
  println("rssi_records = ", rssi_records)
  # println("ssm = $(ssm)")
  println("start estimate path with viterbi")
  @time ep = estimate_path_viterbi(rssi_records, ssm, project.plan )
  # println("end estimate path with viterbi")
  println("path = ", ep)
  push!(estim_path, ep)

  println("Saving for client $(file)")
  save("$(project.path_init_data)/clients_jld/client_$(fileind).jld", "signalrecords", rssi_records,"estimatedpath",ep)
  # save("$(current_path)/res/paths/client_$(fileind).jld", "signalrecords", signal_records,"estimatedpath",ep,"timestamps", time)
end
return estim_path
end


function predprocessing(rssi_p,ap_p,t_p)
  rssi = []
  ap = []
  t = []
  loc_rssi = []
  loc_ap = []
  loc_t = []
    for i=1:length(rssi_p)
      if rssi_p[i]!=0
        push!(rssi,rssi_p[i])
        push!(ap,ap_p[i])
        push!(t,t_p[i])
      end
    end
  # for i=1:length(rssi)
    # println("RSSI = $(rssi[i]) \t ap = $(ap[i]) \t time = $(t[i])")
  # end
  for i=2:length(rssi)
    if t[i]!=t[i-1]
      push!(loc_rssi,rssi[i])
      push!(loc_ap,ap[i])
      push!(loc_t,t[i])
    elseif t[i]==t[i-1] && ap[i]==ap[i-1]
      push!(loc_rssi,(rssi[i]+rssi[i-1])/2)
      push!(loc_ap,ap[i])
      push!(loc_t,t[i])
      i+1
    end
  end
  return loc_rssi, loc_ap, loc_t
end


function get_records(client)
  rssi = Array(client[:,2])
  timestamps = Array{String}(client[:,5])
  ip = Array(client[:,3])
  return rssi,ip,timestamps
end

function getBasicRecords(client)
  rssi = Array(client[:,2])
  timestamps = Array{Int64}(client[:,3])
  ap = Array(client[:,1])
  return rssi,ap,timestamps
end

function convertTimeStamp(timestamps::Array{String})
  time_note = []
  for i=1:length(timestamps)
    push!(time_note,time(Libc.strptime("%Y-%m-%d %T",timestamps[i])))
  end
  println("done convert")
  return time_note
end
function convertIpInID(ip, aps)
  apId = []
  for i=1:length(ip)
      push!(apId,aps[ip[i]])
  end
  return apId
end


end
