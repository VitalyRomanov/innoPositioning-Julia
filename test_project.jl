current_path = pwd()
push!(LOAD_PATH, "$(current_path)/src")

using CoverageMapProject
using JLD
using MapVis

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
else
  println("Unknown choice")
end

# MapVis.visualizeWallVis(proj)

params = [147.55,-20*log10(2.4e9),0.,-0.,-3.,-9.51,-41.14]
# 147.55, -187.604, 0.0, -11.8535, -2.0, -5.0, -15.0 with bounds
# 147.55, -187.604, 0.0, -16.8323, -1.75058, -5.00855, -1.08531 no bounds
params = [147.55, -187.604, 0.0, -11.8535, -2.0, -5.0, -15.0]
if params == []
    params = CoverageMapProject.fit_parameters(proj)
end
CoverageMapProject.calculate_coverage_map(proj,parameters = params,from_dump = false)

ssms = CoverageMapProject.load_ssms(proj)

clients = CoverageMapProject.load_clients(proj)


#warmup
# warmup = CoverageMapProject.load_warmup()
# CoverageMapProject.restore_paths(warmup,ssms,proj)

CoverageMapProject.restore_paths(clients,ssms,proj)

# space = project.plan.limits[1:2,:]
# grid_size = 1.
# grid = convert(Array{Int},floor.((space[:,2] - space[:,1]) / grid_size))
# SpaceInfo(space,grid,grid_size)
