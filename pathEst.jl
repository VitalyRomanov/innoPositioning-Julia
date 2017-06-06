current_path = pwd()
# current_path = "/Users/LTV/Dropbox/mapGen/fast_vis_matr"
# cd(current_path)
push!(LOAD_PATH, "$(current_path)/src")


using HDF5,JLD
using DataFrames
using PyPlot
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
client = readtable("$(current_path)/res/rssi_data/client_data.csv")
# client = readtable("1020.cleaned.csv")
# client = ip_to_labels!(client)
# client = filter_rssi_values!(client,[0,1])
# client = change_rssi_scale!(client)
# client = remove_milliseconds(client)
# client = time_to_ms(client)
# writetable("client_data.csv",client)
# plot_hist_by_apa(client,2)

signal_records = get_rssi_records(client)


# signal_spaces = Array(Array{Float64},1)
# space = [6 6]
# spaceInfo = getSpaceInfo(space,.2)
# pmodel = PModelParam(3.,4.,-53.)
# signal_map = pmodel_space([0 0],spaceInfo,pmodel)
# signal_spaces[1] = signal_map
# test_path,signal_records = path_generation([0.1 0.1],6,signal_map,spaceInfo)
# test_path = [.1 .1; 2. 1.; 4. 3.; 4. 4.; 3. 4.5; 2. 5.]
# signal_records,index_path,index_path_unfolded,vel = path_to_rssi(test_path,signal_map,spaceInfo)


# @time ep,tr = estimate_path_viterbi(signal_records,signal_spaces,spaceInfo,seed = [.1,.1],testing = true, ip = index_path, ipu = index_path_unfolded)

## Estimate path
# test_path = [RssiRecord(-40.,1,0),RssiRecord(-40.,2,1000),RssiRecord(-40.,3,1000),RssiRecord(-40.,4,1000),RssiRecord(-40.,5,1000),
# RssiRecord(-40.,6,1000),RssiRecord(-40.,7,1000),RssiRecord(-40.,8,1000),RssiRecord(-40.,9,1000),RssiRecord(-40.,10,1000)]
# signal_records = test_path
@time ep,tr = estimate_path_viterbi(signal_records,signal_spaces,spaceInfo)


## Display town
walls,bbox = MapPlan.read_data("$(current_path)/res/coverage/walls.txt")
MapPlan.plot_walls(walls)
PyPlot.plot(ep[:,1],ep[:,2],"b")
savefig("$(current_path)/res/recovered_path.png")


for (ind,record) in enumerate(signal_records)
  println("$(record.rssi) $(record.ap) $(record.dt) $(ep[ind,:])")
end

save("ep1.jld", "sr1", signal_records,"ep1",ep)
