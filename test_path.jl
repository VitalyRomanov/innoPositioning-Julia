using Distributions
using Base.Test

struct TestPlan
    limits
    sp2d_size
    gcell_size
end

const rssiSigma = 1 # rssi noise stdev
const sqrt_2 = sqrt(2)

coord2grid(coord,plan) = Int.(ceil.( (coord-plan.limits[1:2,1]) / plan.gcell_size ))
grid2coord(grid_loc,plan) = grid_loc*plan.gcell_size + plan.limits[1:2,1] - plan.gcell_size/2 # .5 is the half of sector size of 1 m2
unfold_index(grid_loc,plan) = grid_loc[1]+(grid_loc[2]-1)*plan.sp2d_size[1]
lognormpdf(x,m,s)  = -((x.-m)./s).^2/2-log(2.5066282746310002.*s)

uniform_location(plan) = begin
   #   creates a distribution over space with all positions beimg equally likely
     spacesize = prod(plan.sp2d_size)
     ones(Float64,tuple(plan.sp2d_size...))/spacesize
 end

 fold_index(index,plan) = begin
      ii = [0,0]
      ii[2] = ceil(index/plan.sp2d_size[1])
      ii[1] = index[1] - (ii[2]-1)* plan.sp2d_size[1]
      ii
  end

function spaceInit(plan)
    g_size = plan.gcell_size # grid resolution
    s_size = plan.sp2d_size
    limits = plan.limits
    global x_grid = collect( limits[1,1]:g_size:(limits[1,2]-1) ) + g_size/2
    global y_grid = collect( limits[2,1]:g_size:(limits[2,2]-1) ) + g_size/2
end

function unit_test_p_rssi(rssi::Float64, # current rssi value
                    signalMap::Array{Float64} # unrolled coverage map
                    )
    # println(typeof(Normal.(signals_Map,rssiSigma)))
    # d = Normal.(signals_Map,rssiSigma)
    # for ind=1:length(d)
        # p_rssi = log.(pdf.(d[ind],rssi))
    # end
    # println("\n\n\n unit_test_p_rssi\n")
    dist = log.(pdf.(Normal.(signalMap,1),rssi))
    # dist = lognormpdf(rssi,signalMap,rssiSigma)
    # println(" dist=\n $(dist)\n\n\n")
    return dist
end
function unit_p_rssi(rssi::Float64, # current rssi value
                    signalMap::Array{Float64} # unrolled coverage map
                    )
    dist = lognormpdf(rssi,signalMap,rssiSigma)
    return dist
end
function run_unnit_tests_p_transition()
    planTest = TestPlan([0 5;0 1;-5 5],[5, 1],1.)
    spaceInit(planTest)
    @test unit_test_p_transition(1,[0. 0.],planTest,1.) == [0.0 0.0; 2.0 0.0; 4.0 0.0; 6.0 0.0; 8.0 0.0]
    @test unit_test_p_transition(2,[2. 0.],planTest,1.) == [-4.0 0.0; -2.0 0.0; 0.0 0.0; 2.0 0.0; 4.0 0.0]
    @test unit_test_p_transition(3,[0. 0.],planTest,1.) == [-4.0 0.0; -2.0 0.0; 0.0 0.0; 2.0 0.0; 4.0 0.0]
    @test unit_test_p_transition(4,[2. 0.],planTest,1.) == [-8.0 0.0; -6.0 0.0; -4.0 0.0; -2.0 0.0; 0.0 0.0]
    @test unit_test_p_transition(5,[0. 0.],planTest,1.) == [-8.0 0.0; -6.0 0.0; -4.0 0.0; -2.0 0.0; 0.0 0.0]

    planTest = TestPlan([0 1;0 5;-5 5],[1, 5],1.)
    spaceInit(planTest)
    @test unit_test_p_transition(1,[0. 0.],planTest,1.) == [0.0 0.0; 0.0 2.0; 0.0 4.0; 0.0 6.0; 0.0 8.0]
    @test unit_test_p_transition(2,[0. 2.],planTest,1.) == [0.0 -4.0; 0.0 -2.0; 0.0 0.0; 0.0 2.0; 0.0 4.0]
    @test unit_test_p_transition(3,[0. 0.],planTest,1.) == [0.0 -4.0; 0.0 -2.0; 0.0 0.0; 0.0 2.0; 0.0 4.0]
    @test unit_test_p_transition(4,[0. 2.],planTest,1.) == [0.0 -8.0; 0.0 -6.0; 0.0 -4.0; 0.0 -2.0; 0.0 0.0]
    @test unit_test_p_transition(5,[0. 0.],planTest,1.) == [0.0 -8.0; 0.0 -6.0; 0.0 -4.0; 0.0 -2.0; 0.0 0.0]

    planTest = TestPlan([0 2;0 2;-5 5],[2, 2],1.)
    spaceInit(planTest)
    @test unit_test_p_transition(1,[0. 0.],planTest,1.) == [0.0 0.0; 2.0 0.0; 0.0 2.0; 2.0 2.0]
    @test unit_test_p_transition(2,[2. 0.],planTest,1.) == [-4.0 0.0; -2.0 0.0; -4.0 2.0; -2.0 2.0]
    @test unit_test_p_transition(3,[0. 2.],planTest,1.) == [0.0 -4.0; 2.0 -4.0; 0.0 -2.0; 2.0 -2.0]
    @test unit_test_p_transition(4,[2. 2.],planTest,1.) == [-4.0 -4.0; -2.0 -4.0; -4.0 -2.0; -2.0 -2.0]

    # @test unit_test_p_transition(6,[0. 0.],planTest,1.) == [0.0 0.0; 2.0 0.0; 4.0 0.0; 6.0 0.0; 8.0 0.0]
    # @test unit_test_p_transition(7,[0. 2.],planTest,1.) == [-4.0 0.0; -2.0 0.0; 0.0 0.0; 2.0 0.0; 4.0 0.0]
    # @test unit_test_p_transition(8,[0. 0.],planTest,1.) == [-4.0 0.0; -2.0 0.0; 0.0 0.0; 2.0 0.0; 4.0 0.0]
    # @test unit_test_p_transition(9,[0. 2.],planTest,1.) == [-8.0 0.0; -6.0 0.0; -4.0 0.0; -2.0 0.0; 0.0 0.0]
    # @test unit_test_p_transition(10,[0. 0.],planTest,1.) ==[-8.0 0.0; -6.0 0.0; -4.0 0.0; -2.0 0.0; 0.0 0.0]
end

function run_unnit_tests_p_rssi(ssm,rssis)
    test_plan = TestPlan([-1 4;0 5;-5 5],[5, 5],1.)
    spaceInit(planTest)
    @test unit_p_rssi(rssis[1].,[0. 0.],planTest,1.) == [0.0 0.0; 2.0 0.0; 4.0 0.0; 6.0 0.0; 8.0 0.0]
    @test unit_p_rssi(2,[2. 0.],planTest,1.) == [-4.0 0.0; -2.0 0.0; 0.0 0.0; 2.0 0.0; 4.0 0.0]
    @test unit_p_rssi(3,[0. 0.],planTest,1.) == [-4.0 0.0; -2.0 0.0; 0.0 0.0; 2.0 0.0; 4.0 0.0]
    @test unit_p_rssi(4,[2. 0.],planTest,1.) == [-8.0 0.0; -6.0 0.0; -4.0 0.0; -2.0 0.0; 0.0 0.0]
    @test unit_p_rssi(5,[0. 0.],planTest,1.) == [-8.0 0.0; -6.0 0.0; -4.0 0.0; -2.0 0.0; 0.0 0.0]



function unit_test_p_transition(state,
                                v::Array{Float64}, # velocity defines mean value
                                plan_test,
                                dt)
    # println("\n\n\n unit_test_p_transition \n\n\n")
    g_size = plan_test.gcell_size # grid resolution
    s_size = plan_test.sp2d_size
    limits = plan_test.limits
    loc = grid2coord(fold_index(state,plan_test),plan_test)
    dx = x_grid #  * g_size + g_size/2
    dy = y_grid #  * g_size + g_size/2
    # x = normalpdf(dx,loc[1]+v[1]*dt,dt/sqrt_2)
    # y = normalpdf(dy,loc[2]+v[2]*dt,dt/sqrt_2)
    x = pdf.(Normal.(loc[1]+v[1]*dt,dt/sqrt_2),dx)
    y = pdf.(Normal.(loc[2]+v[2]*dt,dt/sqrt_2),dy)
    brc = broadcast(*,x,y')
    if sum(brc) == 1
        println("sum == 1")
    end
    brc ./= sum(brc)
    brc .= log.(brc)
    ax = (dx-loc[1]-v[1]*dt)*2/dt^2
    ay = (dy-loc[2]-v[2]*dt)*2/dt^2

    vx = (ones(Float64,s_size[2])*(v[1]+ax*dt)')'
    vy = ones(Float64,s_size[1])*(v[2]+ay*dt)'
    # vx = (ones(Float64,s_size[2])*((2*dx/dt)-v[1])')'
    # vy = ones(Float64,s_size[1])*((2*dy/dt)-v[2])'
    # return brc,[vx[:] vy[:]]
    return [vx[:] vy[:]]
end



function unit_test_cp(p_rssi,step)
    # println("\n\n\n unit_test_cp \n")
    # println("p_rssi=\n$(p_rssi)\n\n")
  # println("step=\n$(step)\n\n")
  # rssi = exp.(p_rssi)
  # steps = exp.(step)
  # step ./= sum(step)
  # prob = rssi + steps
  prob = p_rssi+step
  prob.= exp.(prob)
  prob./=sum(prob)
  prob.=log.(prob)
  # prob .= log.(prob)
  # println("prob=$(prob)")
  return prob[:] # should do normalization?
end

function unit_test_ct(trellis_part,velocity_part,p_transition,v_transition,path_back,state)
    # println("\n\n\n unit test_ct \n\n\n")
  for ss in 1:1:length(trellis_part)
      # println("\n ss=$(ss) \n")

      # println("\n p_transition=$(p_transition[ss]) \n")
      # println("\n trellis_part=$(trellis_part[ss]) \n")
    if p_transition[ss] > trellis_part[ss]# maybe >=
      trellis_part[ss] = p_transition[ss]
      velocity_part[ss,:] = v_transition[ss,:]
      # println("\n state=$(state) \n")
      path_back[ss] = state
    elseif p_transition[ss] == trellis_part[ss]
        println("state = $(state)")
        println("ss = $(ss)")

      println("p_transition = $(p_transition[ss]) and Trellis = $(trellis_part[ss])\n\n\n")

    end
  end
  return trellis_part,velocity_part,path_back
end

function unit_beam_search(data, K) # data = trellis[:,1]
  # println("\n\n\n unit_beam_search \n\n\n")
  trellis = sort(data, rev=true)
  state = findin(data,trellis[1:K])
  # println("state from beam = $(state)")
  return state
end

function unit_test_viterbi(signals,
                                signal_map,
                                plan;
                                seed = [-1. -1.])
  # plan,signal_map = o_plan,o_signal_map
  spaceInit(plan)
  g_size = prod(plan.sp2d_size) # total number of cells in the grid
  trellis = zeros(Float64,g_size,2)-1e10 # stores previous and current state probability
  velocity = zeros(Float64,g_size,2,2) # stores previous and current state velocity in xy dimensions
  path_back = zeros(Int32,g_size,length(signals)) # stores feasible paths
  init_step = zeros(Float64,tuple(plan.sp2d_size...))
  if seed==[-1. -1.]
    init_step = uniform_location(plan)
  else
    init_step = seed_location(seed,plan)#*1e100
  end
  c_rssi = signals[1].rssi
  c_ap_id = signals[1].ap
  p_rssi = unit_test_p_rssi(c_rssi,
                        signal_map[c_ap_id]) # the value is boosted by exp(100)
  trellis[:,1] = unit_test_cp(p_rssi,init_step) # values sum to 1e100?
  println("trellis1 $(1) = \n $(trellis[:,1])\n\n ")
  println("prob $(1) = \n $(findmax(trellis[:,1]))\n\n ")

  # trellis[:,1] .= exp.(trellis[:,1])
  # trellis[:,1] ./= sum(trellis[:,1])
  # trellis[:,1] .= log.(trellis[:,1])
  # println("trellis1 $(1) = \n $(trellis[:,1])\n\n ")
  # println("prob $(1) = \n $(findmax(trellis[:,1]))\n\n ")

  for i = 1:1:length(signals)-1
    dt = signals[i+1].t - signals[1].t
    if dt==0.
      for ii=1:length(path_back[:,i+1])
        path_back[ii,i+1] = ii
      end
      continue
    end

    c_rssi = signals[i+1].rssi
    c_ap_id = signals[i+1].ap
    p_rssi = unit_test_p_rssi(c_rssi,
                          signal_map[c_ap_id])
    # for state in 1:1:g_size
    for state in unit_beam_search(trellis[:,1],25)
      # if (trellis[state,1]>-100)
        p_transition,v_transition = unit_test_p_transition(state,
                                                velocity[state,:,1],
                                                plan,
                                                dt)
        cp = unit_test_cp(p_rssi,p_transition)

        trellis[:,2],velocity[:,:,2],path_back[:,i+1] = unit_test_ct(trellis[:,2],
                                                                        velocity[:,:,2],
                                                                        cp+trellis[state,1],
                                                                        v_transition,
                                                                        path_back[:,i+1],
                                                                        state)

     if (i+1) == 7

         println("v_transition $(i+1) state = $(state) = \n $(v_transition)\n\n ")
     end
      # end
    end
    println("trellis2 $(i+1) = \n $(trellis[:,2])\n\n ")
    println("prob $(i+1) = \n $(findmax(trellis[:,2]))\n\n ")
    println("path_back $(i+1) = \n $(path_back[:,i+1])\n\n ")
    println("velocity $(i+1) = \n $(velocity[:,:,2])\n\n ")
    # trellis[:,2] .= exp.(trellis[:,2]) -> X
    # trellis[:,2] ./= sum(trellis[:,2]) -> X
    # trellis[:,2] .= log.(trellis[:,2]) -> X
    # println("\n Trellis 2 after notmalization\n\n ")
    # println("trellis2 $(i+1) = \n $(trellis[:,2])\n\n ")
    # println("prob $(i+1) = \n $(findmax(trellis[:,2]))\n\n ")
    # lpdfmax = maximum(trellis[:,2])
    # println("\n\n\nThe maximum value of trellis = $(lpdfmax) in o =$(i+1)")
    # println("\n\n\nindex of The maximum value trellis=\n\n $(find(a->a==lpdfmax, trellis[:,2]))\n\n\n\n")
    trellis[:,1] = trellis[:,2]
    velocity[:,:,1] = velocity[:,:,2]
    trellis[:,2] = -1e10 #zeros(Float64,g_size)-1e10 # creating new array -> performance issues
    velocity[:,:,2] = 0. #zeros(Float64,g_size,2) # creating new array -> performance issues
  end
  estimated_path = zeros(Float64,length(signals),3)
  step_back = indmax(trellis[:,1])
  # println(trellis[:,1])
  # println(findmin(trellis[:,1]))

  # println(step_back)
  steps = []
  for step_ind in length(signals):-1:1
    push!(steps,step_back)
    estimated_path[step_ind,1:2] = grid2coord(
                                    fold_index(step_back,plan),
                                                plan)
    estimated_path[step_ind,3] = signals[step_ind].t
    step_back = path_back[step_back,step_ind]


  end
  return estimated_path, steps , path_back
end



ap = [1.2, 1.8]
path = [.6 .6;
 1.9 .6;
 2.4 .7;
 3.1 .9;
 3.6 2.3;
 3.5 3.6;
 2.7 3.2;
 1.8 2.9;
 1.7 2.8;
 1.3 2.6;
 .6 2.2]

 t = [ 0.0, 1.5, 2.0, 2.9, 4.5, 5.7, 7.8, 9.4, 9.6, 10.0, 10.8]


n = size(path,1)

mutable struct RssiRecord
    rssi
    ap
    t
end

rssis = Array{RssiRecord}(n)


rssi(d) = begin
    d0 = 1
    p_d0 = -10
    alpha = 3.1
    val = p_d0 - 10 * alpha * log10(d/d0)
    val
end

for i=1:n
    d = norm(ap-path[i,:])
    d0 = 1
    p_d0 = -10
    alpha = 3.1
    val = rssi(d)
    rr = RssiRecord(val,1,t[i])
    rssis[i] = rr
    println("t$(i) = $(val)")
end

ssm = Array{Float64}(5,5)

for (x_i,x) = enumerate(-.5:1:3.55), (y_i,y) = enumerate(.5:1:5)
    d = norm(ap-[x,y])
    ssm[x_i,y_i] = rssi(d)
    println("ssm x=$x y=$y rssi=$(ssm[x_i,y_i])")
end

test_plan = TestPlan([-1 4;0 5;-5 5],[5, 5],1.)

planTest = TestPlan([0 5;0 1;-5 5],[5, 1],1.)
ssmTest = Array{Float64}(5,1)

for (x_i,x) = enumerate(0.5:1:5), (y_i,y) = enumerate(0:1:0)
    d = norm(ap-[x,y])
    ssmTest[x_i,y_i] = rssi(d)
    println("ssm x=$x y=$y rssi=$(ssm[x_i,y_i])")
end

# spaceInit(test_plan)

# c_rssi = rssis[1].rssi
# p_rssi = unit_test_p_rssi(c_rssi,
                      # ssm) # the value is boosted by exp(100)




ep_path,ep_state,back = unit_test_viterbi(rssis,[ssm],test_plan)




d = sqrt((1.2-0.6)^2 + (1.8-0.6)^2)
rssi(d)

using JLD
JLD.save("test_ssm_1.jld","ssm",ssm)
JLD.save("test_path.jld","path",rssis,"gt",path)
