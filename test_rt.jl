current_path = "/Users/LTV/dev_projects/innoPositioning-Julia"
push!(LOAD_PATH, "$(current_path)/src")

using MapPrimitives
using MapPlan
using RadixTree

function wall2mbr(wall::Wall)
  v1 = minimum([wall.points[1].val wall.points[2].val wall.points[3].val wall.points[4].val],2)
  v2 = maximum([wall.points[1].val wall.points[2].val wall.points[3].val wall.points[4].val],2)
  return MBR(v1,v2)
end

walls,x_lims,y_lims = MapPlan.read_data("$(current_path)/res/coverage/walls_tables.txt")

lims = [x_lims';y_lims';[0,6]']

tree = RadixTree.create_radix_tree(RadixTree.obj2mbr(walls,wall2mbr),lims)

overlap = RadixTree.probe(tree,RadixTree.MBR([30.0006,10.0008,1.0],[34.9994,15.9992,1.0]))

println(overlap)
