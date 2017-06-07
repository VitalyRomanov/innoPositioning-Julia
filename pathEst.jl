current_path = pwd()
# current_path = "/Users/LTV/Dropbox/mapGen/fast_vis_matr"
# cd(current_path)
push!(LOAD_PATH, "$(current_path)/src")
data_folder = "/Users/LTV/Documents/rssi_data/client_data/"


using HDF5,JLD
using DataFrames
using ClientDataProc
using MapPrimitives
using LocTrack
using MapPlan


# read signal distribution matrices
signal_spaces = Array{Float64}[load("$(current_path)/res/coverage/data_$(i).jld","data") for i in 1:10]

print("Coverage for $(length(signal_spaces)) is loaded\n")

# Set up map space
spaceInfo = getSpaceInfo([size(signal_spaces[1],1),size(signal_spaces[1],2)],2.)

signal_spaces = downsample_maps(signal_spaces,2)

# Read and filter
# ClientDataProc.filter_and_process("/Users/LTV/Documents/rssi_data/1020.cleaned.csv")

# client = readtable("/Users/LTV/Documents/rssi_data/1020.cleaned.csv")
# client = ip_to_labels!(client)
# client = filter_rssi_values!(client,[0,1])
# client = change_rssi_scale!(client)
# client = remove_milliseconds(client)
# client[:,5] = deepcopy(client[:,4])
# client = ClientDataProc.time_to_ms(client)
# writetable("/Users/LTV/Documents/rssi_data/client_data/client_01.csv",client)
# plot_hist_by_apa(client,2)



for (fileind,file) in enumerate(readdir(data_folder))
  println("Reading file $(file)")
  client = readtable("$(data_folder)/$(file)")
  println("Formatting...")
  signal_records,timestamps = get_rssi_records(client)
  ep,~ = estimate_path_viterbi(signal_records,signal_spaces,spaceInfo)
  println("Saving")
  save("$(current_path)/res/paths/client_$(fileind).jld", "signalrecords", signal_records,"estimatedpath",ep,"timestamps", timestamps)
end

timestamps = Array(client[:,5])
timestamps = timestamps[1:11922]
save("$(current_path)/res/paths/client_1.jld", "signalrecords", signal_records,"estimatedpath",ep,"timestamps", timestamps)

# client = readtable("/Users/LTV/Documents/rssi_data/client_data/client_01.csv")

# signal_records = get_rssi_records(client)


# @time ep,tr = estimate_path_viterbi(signal_records,signal_spaces,spaceInfo)


## Display town
# walls,bbox = MapPlan.read_data("$(current_path)/res/coverage/walls.txt")
# MapPlan.plot_walls(walls)
# PyPlot.plot(ep[:,1],ep[:,2],"b")
# savefig("$(current_path)/res/recovered_path.png")


# for (ind,record) in enumerate(signal_records)
#   println("$(record.rssi) $(record.ap) $(record.dt) $(ep[ind,:])")
# end

# save("$(current_path)/res/paths/client_01.jld", "sr1", signal_records,"ep1",ep,"timestamp", Array(client[:,5]))
