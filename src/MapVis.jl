module MapVis

using Plots
using StatPlots
using StatsBase

gr()

function plot_walls!(project; ids = [], col = :black, lw = 2)
  if ids == []
      ids = 1:length(project.plan.walls)
  end
  for wall in project.plan.walls[ids]
    if wall.polygon[1][3]==wall.polygon[3][3]
      continue
    end
    xs = [wall.polygon[1][1],wall.polygon[4][1]]-project.plan.limits[1,1]
    ys = [wall.polygon[1][2],wall.polygon[4][2]]-project.plan.limits[2,1]
    plot!(xs,ys,
        linecolor=col,
        linewidth=lw,
        # xlims = project.plan.limits[1,:],
        # ylims = project.plan.limits[2,:],
        legend = false)
  end
end



function plot_map(ssm,project,map_ind,filename)
    rssi_min = -200.
    rssi_max = maximum(ssm)

    println("Minimum: $(rssi_min)   Maximum: $(rssi_max)")

    plot(ssm,
        seriestype=:heatmap,
        seriescolor=ColorGradient([colorant"white", colorant"orange", colorant"red"]),
        zlims=(rssi_min,rssi_max),
        legend = false,
        grid=false,
        axis=false)

    plot_walls!(project)

    # Plot the marker for AP location
    plot!([project.APs[map_ind][1]-project.plan.limits[1,1]],
        [project.APs[map_ind][2]-project.plan.limits[2,1]],
        markershape=:diamond,
        markercolor=:pink)

    savefig(filename)
    # savefig("map_.png")
end


function plot_map(project,map_ind)
  ssm = project.ssms[map_ind]'
  ssm_map_path = "$(project.path_init_data)/map_$(map_ind).png"

  plot_map(ssm,project,map_ind,filename)
end


function plot_paths(ssm, project, real_path, est_path, point = false, real = false, ind = 1)
  ssm_map_path = "$(project.path_init_data)/map_$(ind)"

  plot_map(ssm, project,length(project.APs),ssm_map_path)

  limx = project.plan.limits[1,1]
  limy = project.plan.limits[2,1]
  println("limx = $(limx) , limy = $(limy)")
  for path in real_path
    for i=1:size(path,1)-1
      plot!([path[i,1],path[i+1,1]]-limx,[path[i,2],path[i+1,2]]-limy,linecolor=:green, subplot = 1)
    end
  end
  for path in est_path
    for i=1:size(path,1)-1
      plot!([path[i,1],path[i+1,1]]-limx,[path[i,2],path[i+1,2]]-limy,linecolor=:red, subplot = 1)
    end
  end
  if point
    for path in est_path
      for i=1:size(path,1)
        plot!([path[i,1]-limx],[path[i,2]-limy],markershape=:diamond,markercolor=:red,subplot = 1)
        annotate!((path[i,1]-limx),(path[i,2]-limy)+2*7/10,text("$i",7),subplot = 1)
      end
    end
    for path in real_path
      for i=1:size(path,1)
        plot!([path[i,1]-limx],[path[i,2]-limy],markershape=:diamond,markercolor=:green,subplot = 1)
        annotate!((path[i,1]-limx),(path[i,2]-limy)+2*7/10,text("$i",7),subplot = 1)
      end
    end
  end
  for map_ind=1:length(project.APs)
      plot!([project.APs[map_ind][1]-limx], [project.APs[map_ind][2]-limy], markershape=:diamond, markercolor=:pink)
      annotate!((project.APs[map_ind][1]-limx),(project.APs[map_ind][2]-limy)+2*7/10,text("$(map_ind)",7),subplot = 1)
  end
    if !isdir("$(project.path_init_data)/paths")
      mkdir("$(project.path_init_data)/paths")
    end
    savefig("estimated_path_$(ind).svg")
    savefig("$(project.path_init_data)/paths/estimated_path_$(ind).svg")
  end

  function plot_paths(ssm, project, est_path, ind)
    ssm_map_path = "$(project.path_init_data)/map_$(ind)"

    plot_map(ssm, project,length(project.APs),ssm_map_path)
    limx = project.plan.limits[1,1]
    limy = project.plan.limits[2,1]
    println("limx = $(limx) , limy = $(limy)")
    for path in est_path
      for i=1:size(path,1)-1
        plot!([path[i,1],path[i+1,1]]-limx,[path[i,2],path[i+1,2]]-limy,linecolor=:red, subplot = 1)
      end
    end
      for path in est_path
        for i=1:size(path,1)
          plot!([path[i,1]-limx],[path[i,2]-limy],markershape=:diamond,markercolor=:red,subplot = 1)
          annotate!((path[i,1]-limx),(path[i,2]-limy)+2*7/10,text("$i",7),subplot = 1)
        end
      end
    for map_ind=1:length(project.APs)
        plot!([project.APs[map_ind][1]-limx], [project.APs[map_ind][2]-limy], markershape=:diamond, markercolor=:pink)
        annotate!((project.APs[map_ind][1]-limx),(project.APs[map_ind][2]-limy)+2*7/10,text("$(map_ind)",7),subplot = 1)
    end
      if !isdir("$(project.path_init_data)/paths")
        mkdir("$(project.path_init_data)/paths")
      end
      savefig("$(project.path_init_data)/paths/estimated_path.svg")
    end


function plot_paths(project,paths,mp_ind)
  plot(size(600,600),
      axis = false)
  plot_walls!(project)

  limx = project.plan.limits[1,1]
  limy = project.plan.limits[2,1]

  for path in paths
    for i=1:length(path)-1
      plot!([path[i][1],path[i+1][1]]-limx,[path[i][2],path[i+1][2]]-limy,linecolor=:red)
    end
  end

  if !isdir("$(project.path_init_data)/paths")
    mkdir("$(project.path_init_data)/paths")
  end

  savefig("$(project.path_init_data)/paths/$(mp_ind).svg")
end



function visualizeWallVis(project, range = [])
    # if range not specified plot for all walls
    if range == []
        range = 1:length(project.plan.walls)
    end
    println("Visualizing walls visibility in range ", range)
    # Path for storing visibility visualizations
    wallVisOutPath = "$(project.path_init_data)/wall_Vis"
    # Create a folder to store visualizations
    if !isdir(wallVisOutPath)
        mkdir(wallVisOutPath)
    end
    # Create figure for every wall in range
    for wall_id = range
        print("\rVisualizing visibility for wall $(wall_id).........")
        plot()
        plot_walls!(project)
        plot_walls!(project,
                  ids = find(project.plan.vis_matr[wall_id,:]),
                  col = :red,
                  lw=4)
        plot_walls!(project,
                  ids = [wall_id],
                  col = :green,
                  lw=4)
        savefig("$(wallVisOutPath)/wall_$(wall_id).svg")
    end
    println("done")
end


function plot_fit_error(project,dist,mean_rssi,ap_id,measur)
    # create a plot with specified parameters
    # this plot will contain the area map in the first subplot and
    # comparison of the observed and estimated rssi values in the second
    fs = 7 # font size for plot annotations
    plot(legend = false,
        grid=false,
        # axis=false,
        size=(900,1800),
        layout=grid(2,1,heights=[.5,.5]))
    plot_walls!(project)


    order = sortperm(dist)

    plot!(dist[order],mean_rssi[order],legend=false,subplot=2)
    plot!([project.APs[ap_id][1]-project.plan.limits[1,1]],[project.APs[ap_id][2]-project.plan.limits[2,1]],
        markershape=:diamond,
        markercolor=:pink,
        subplot = 1)

    for elem in order
      x = measur[elem].location[1]-project.plan.limits[1,1]
      y = measur[elem].location[2]-project.plan.limits[2,1]
      plot!([x],[y],
          markershape=:diamond,
          markercolor=:blue,
          subplot = 1)
      annotate!(x,y+2*fs/10,text("$elem",fs),subplot = 1)
      boxplot!([dist[elem]],measur[elem].rssi,lab="Loc $(elem)",subplot = 2)
      annotate!(dist[elem]*1.01,mean(measur[elem].rssi)*1.01,text("$elem"),subplot = 2)
    end

    savefig("$(project.path_init_data)/all_it_takes_$(ap_id).svg")
end

end
