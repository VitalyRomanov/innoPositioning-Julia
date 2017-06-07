current_path = pwd()
# current_path = "/Users/LTV/dev_projects/innoPositioning-Julia"
push!(LOAD_PATH, "$(current_path)/src")

import HDF5,JLD
using LocTrack
using MapPlan
using Plots


~,bbox = MapPlan.read_data("$(current_path)/res/coverage/walls.txt")
x_max = bbox[3].points[1].val[1]
y_max = bbox[3].points[1].val[2]

walls = MapPlan.read_data_2d("$(current_path)/res/coverage/walls2d.txt")

println("Drawing the town")
# town_plan = MapPlan.plot_walls_2d(walls,x_max,y_max)

data_folder = "$(current_path)/res/paths"

for (fileind,file) in enumerate(readdir(data_folder))
  println("Reading $(file)")
  user_dump = JLD.load("$(data_folder)/$(file)")
  ep = user_dump["estimatedpath"]
  signal_records = user_dump["signalrecords"]
  timestamps = user_dump["timestamps"]

  nn = 0.
  cw = 0.

  gr()

  anim = @animate for i=2:100length(timestamps)
    MapPlan.plot_walls_2d(walls,x_max,y_max)
    plot!([ep[i-1,1],ep[i,1]],[ep[i-1,2],ep[i,2]],linecolor = :blue,xlims=(0,x_max),ylims=(0,y_max),grid=false,legend=false,axis=false)
    plot!([ep[i,1]],[ep[i,2]],markershape=:circle,markercolor=:red,xlims=(0,x_max),ylims=(0,y_max),annotations=(x_max-420,20,text("$(timestamps[i])",:left,:red,16,"Courier Bold")),grid=false,legend=false,axis=false)

    nn = i/length(timestamps)*100
    if nn-cw>1.
      print("\r$(i/length(timestamps)*100)           ")
      cw = nn
    end
  end
  gif(anim, "$(data_folder)/$(file[1:end-3]).gif", fps = 30)
end


# user_dump = JLD.load("$(data_folder)/client_1.jld")
# ep = user_dump["estimatedpath"]
# signal_records = user_dump["signalrecords"]
# timestamps = user_dump["timestamps"]
#
#
#
# anim = @animate for i=2:10
#   MapPlan.plot_walls_2d(walls,x_max,y_max)
#   plot!([ep[i-1,1],ep[i,1]],[ep[i-1,2],ep[i,2]],linecolor = :blue,xlims=(0,x_max),ylims=(0,y_max),grid=false,legend=false,axis=false)
#   plot!([ep[i,1]],[ep[i,2]],markershape=:circle,markercolor=:red,xlims=(0,x_max),ylims=(0,y_max),annotations=(x_max-420,20,text("$(timestamps[i])",:left,:red,16,"Courier Bold")),grid=false,legend=false,axis=false)
#   print("\r$(i)")
# end
# gif(anim, "$(data_folder)/$(file[1:end-3]).gif", fps = 30)




# using Plots
# pyplot()
# anim = Animation()
