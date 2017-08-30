current_path = pwd()#"/Users/LTV/dev_projects/innoPositioning-Julia"
# current_path = "/home/vromanov/dev/innoPositioning-Julia"
# current_path = "/Users/LTV/dev_projects/innoPositioning-Julia"
cd(current_path)
push!(LOAD_PATH, "$(current_path)/src")

using CoverageMapProject
using Plots
using JLD


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
  # resp = 3
  resp = parse(Int,readline())
end


if resp==1
  print("Enter new project name: ")
  name = readline()[1:end-1]
  print("Enter path for initial data: ")
  load_path = strip(readline())
  print("Enter saving location: ")
  save_path = strip(readline())
  # name = "town2"
  # load_path = "/home/ltv/Dropbox/work/coverage/town"
  # save_path = "/home/ltv/Dropbox/work/coverage/town"
  proj = CoverageMapProject.create_project(load_path,save_path,name)
  CoverageMapProject.save_session("$(save_path)/$(name).jld")
  # CoverageMapProject.calculate_image_trees(proj)
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



# proj.image_trees_ready = false
# CoverageMapProject.calculate_image_trees(proj)

# CoverageMapProject.calculate_coverage_map(proj)
# CoverageMapProject.plot_map(proj,1)

# params = CoverageMapProject.fit_parameters(proj,1)
# JLD.save("/home/ltv/Documents/coverage/town/vm.jld","vm",proj.plan.vis_matr)
# proj.ssms_ready_count = 0
params = [147.55,-20*log10(2.4e9),0.,-0.,-2.5,-12.53,-100.]
CoverageMapProject.calculate_coverage_map(proj,parameters = params)
# CoverageMapProject.recalculate_coverage_map(proj,1,parameters = params)
# CoverageMapProject.plot_map(proj,1)

# plot(proj.ssms[1]',seriestype=:heatmap,seriescolor=ColorGradient([colorant"white", colorant"orange", colorant"red"]),zlims=(-40,30),legend = false,grid=false,axis=false)
