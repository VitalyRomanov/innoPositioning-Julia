current_path = "/home/ltv/dev/innoPositioning-Julia"
cd(current_path)
push!(LOAD_PATH, "$(current_path)/src")


using Geometry
using MapPrimitives

ray1 = Line([0.,0.,0.],[10.,10.,10])

awall = Wall3D(1,[[0.,0.,0.],[0.,10,0],[10.,10,0],[10.,0,0]])

@time line2mbr(ray1)
