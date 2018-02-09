module LocTrack

  # dt = 1
  using MapPlan
  using Base.Test
  using Distributions
  const rssiSigma = 2 # rssi noise stdev

  export getSpaceInfo,estimate_path_viterbi,RssiRecord,pmodel_space,path_to_rssi,path_generation

  type RssiRecord
    rssi::Float64
    ap::Int
    t::Float64
  end

  # type SpaceInfo
  #   space::Array{Int}
  #   grid::Array{Int}
  #   grid_size::Float64
  # end

  # type PModelParam
  #   alpha::Float32
  #   d0::Float32
  #   p_d0::Float32
  # end

  struct TestPlan
      limits
      sp2d_size
  end

  ##### Auxiliary functions

  # log of Normal PDF
  # lognormpdf(x,m,s)  = @. -((x-m)/s)^2/2-log(2.5066282746310002*s)
  lognormpdf(x,m,s)  = -((x.-m)./s).^2/2-log(2.5066282746310002.*s)
  # lognormpdf(x,m,s)  = @. -((x-m)/s)^2/2-log(2.5066282746310002*s)


  normalpdf(x,m,s)= pdf.(Normal(m,s),x))
  normalpds(x,m,s) = (1/(2.506628274631001.*s)).^(((x.-m)./s).^2/2)
  # function getSpaceInfo(space,grid_size=1.)
  #   grid = Int.(floor.(space/grid_size))
  #   return SpaceInfo(space,grid,grid_size)
  # end

  # function pmodel_space(ap_loc,spaceInfo::SpaceInfo,pmodel::PModelParam)
  #   # creates coverage map
  #   x = ones(Float64,spaceInfo.grid[1],1) * collect(0:1:(spaceInfo.grid[2]-1))' * spaceInfo.grid_size
  #   y = (ones(Float64,spaceInfo.grid[2],1) * collect(0:1:(spaceInfo.grid[1]-1))')' * spaceInfo.grid_size
  #   d = sqrt((x-ap_loc[1]).^2+(y-ap_loc[2]).^2)
  #   d[1,1] = .1
  #   value = pmodel.p_d0 - 10 * pmodel.alpha * log10(d/pmodel.d0)
  #   return value
  # end

  # coord2grid takes location on the map and converts it to grid coords
  coord2grid(coord,plan) = Int.(ceil.(coord-plan.limits[1:2,1]))

  # function coord2grid(coord::Array{Float64},grid_size::Float64)
  #   # covert coordinates on the map to grid index
  #   return Int.(ceil.(coord/grid_size))
  # end

  # grid2coord takes grid coords and converts it to location on the map
  grid2coord(grid_loc,plan) = grid_loc+plan.limits[1:2,1]-.5 # .5 is the half of sector size of 1 m2

  # function grid2coord(index::Array{Int},grid_size::Float64)
  #   # coordinates on map, based on 2d grid index
  #   return index*grid_size-grid_size/2
  # end

  # unfold_index taked 2d grid location and does mapping from 2d to 1d index
  unfold_index(grid_loc,plan) = grid_loc[1]+(grid_loc[2]-1)*plan.sp2d_size[1]

  # function unfold_index(index::Array{Int},spaceInfo::SpaceInfo)
  #   return index[1]+(index[2]-1)*(spaceInfo.grid[1])
  # end

  # fold_index taked 1d index and does mapping from 1d to 2d grid location
  fold_index(index,plan) = begin
      ii = [0,0]
      ii[2] = ceil(index/plan.sp2d_size[1])
      ii[1] = index[1] - (ii[2]-1)* plan.sp2d_size[1]
      ii
  end

  # function fold_index(index,spaceInfo)
  #   ii = zeros(Int64,2,1)
  #   ii[2] = ceil(index/spaceInfo.grid[1])
  #   ii[1] = index[1] - (ii[2]-1)* spaceInfo.grid[1]
  #   return ii
  # end

  seed_location(seed,plan) = begin
    # creates distrubution over space with all probability mass in the position
    # specified by seed
      index = coord2grid(seed,plan)
      seeded_location = zeros(Float64,tuple(plan.sp2d_size...))
      seeded_location[index[1],index[2]] = 1.
      seeded_location
  end

  # function seed_location(seed,spaceInfo)
  #   index = coord2grid(seed,spaceInfo.grid_size)
  #   seeded_location = zeros(Float64,tuple(spaceInfo.grid...))
  #   # seeded_location = zeros(Float64,spaceInfo.grid[1],spaceInfo.grid[2])
  #   seeded_location[index[1],index[2]] = 1.
  #   return seeded_location
  # end

  uniform_location(plan) = begin
    #   creates a distribution over space with all positions beimg equally likely
      spacesize = prod(plan.sp2d_size)
      ones(Float64,tuple(plan.sp2d_size...))/spacesize
  end
  # function uniform_location(spaceInfo)
  #   spacesize = prod(spaceInfo.grid)
  #   return ones(Float64,tuple(spaceInfo.grid...))/spacesize
  #   # return ones(Float64,spaceInfo.grid[1],spaceInfo.grid[2])/spacesize
  # end

  function run_tests()
      plan = TestPlan([-5 5;-5 5;-5 5],[10, 10])

      @test coord2grid([4.9,4.9],plan) == [10,10]
      @test grid2coord([1,1],plan) == [-4.5,-4.5]
      @test unfold_index([1,1],plan) == 1
      @test unfold_index([1,2],plan) == 11
      @test fold_index(1,plan) == [1,1]
      @test fold_index(11,plan) == [1,2]
  end


  ####### Project functions
  # type SpaceInfo
  #  -    space::Array{Int}
  #  -    grid::Array{Int}
  #  -    grid_size::Float64
  #  -  end
 #  function getSpaceInfo(space,grid_size=1.)
 # +  #   grid = Int.(floor.(space/grid_size))
 # +  #   return SpaceInfo(space,grid,grid_size)
 #    # end
 # signal_map::Array{Array{Float64}},
 # plan;
 # seed = [-1. -1.],
 # grid_space_size = spaceInfo.grid[1]*spaceInfo.grid[2] # total number of cells in the grid
 # +    g_size = prod(plan.sp2d_size) # total number of cells in the grid
# grid_size=1.
 function path_generation(distribution, path_lenght,
                          signal_map,
                          plan,
                          seed = [-1. -1.])
  # println("Start path generation")
  path, signals = distribution(path_lenght, signal_map, plan, seed)
  # println("Finish path generation")

  return path,signals
end

function acelerationDist(path_lenght,
                         signal_map,
                         plan,
                         seed,
                         dt = 1.)
  r_sigma = 2 #instead 0.1
  # seed : initial location
  # path_lenght : --||--
  # signal_map : signal strength matrix
  # dt :
  # println("Start acelerationDist")
  t = dt
  path = zeros(Float64,path_lenght,2) # stores generated path here
  # rssi = zeros(Float64,path_lenght,1)
  signals  = Array(RssiRecord,path_lenght)
  v = zeros(Float64,2,1)
  path[1,:] = seed

  index = coord2grid(path[1,:],plan)
  # print(index)
  # rssi[1] = randn()+signal_map[index[1],index[2]]
  rssi = randn()+signal_map[1][index[1],index[2]]
  signals[1] = RssiRecord(rssi,1,1.)
  # println("Start loop to path_lenght")
  for i in 2:1:path_lenght
    valid_sample = false
    s = seed
    a = seed
    t+=1
    # println("Start loop !valid_sample")
    a_stdev = 0.5487618458
    while !valid_sample
      a = randn(2,1) * a_stdev
      s = path[i-1,:] + v*dt + a*dt^2/2
      # println("a =$(a) & v=$(v) & s = $(s)")
      # if (s[1]>0) && (s[2]>0) && (s[1]<spaceInfo.space[1,1]) && (s[2]<spaceInfo.space[2,1])
      # if prod(s.>0) && prod(s.>plan.limits[1:2,1]) && prod(s.<plan.limits[1:2,2])
      if prod(s.>plan.limits[1:2,1]) && prod(s.<plan.limits[1:2,2])
        valid_sample = true
      else
        # println(v, s)
        a_stdev += .05
      end
    end
    # print("Finish loop valid simple")
    path[i,:] = s
    v += a*dt
    index = coord2grid(path[i,:],plan)
    rssi = r_sigma*randn()+signal_map[1][index[1],index[2]]
    signals[i] = RssiRecord(rssi,1,t)
  end
  return path,signals
end

  #  function path_generation(seed,path_lenght,signal_map,spaceInfo, dt = 1.)
  #   # seed : initial location
  #   # path_lenght : --||--
  #   # signal_map : signal strength matrix
  #   # spaceInfo :
  #   # dt :
  #   path = zeros(Float64,path_lenght,2) # stores generated path here
  #   # rssi = zeros(Float64,path_lenght,1)
  #   signals  = Array(RssiRecord,path_lenght)
  #   v = zeros(Float64,2,1)
  #   path[1,:] = seed
  #
  #   index = coord2grid(seed,spaceInfo.grid_size)
  #   # print(index)
  #   # rssi[1] = randn()+signal_map[index[1],index[2]]
  #   rssi = randn()+signal_map[index[1],index[2]]
  #   signals[1] = RssiRecord(rssi,1,1.)
  #
  #   for i in 2:1:path_lenght
  #     valid_sample = false
  #     s = seed
  #     a = seed
  #     while !valid_sample
  #       a = randn(2,1)
  #       s = path[i-1,:] + v*dt + a*dt^2/2
  #       # if (s[1]>0) && (s[2]>0) && (s[1]<spaceInfo.space[1,1]) && (s[2]<spaceInfo.space[2,1])
  #       if prod(s.>0) && prod(s.<spaceInfo.space[1:2,1])
  #         valid_sample = true
  #       end
  #     end
  #     path[i,:] = s
  #     v += a*dt
  #     index = coord2grid(path[i,:],spaceInfo.grid_size)
  #     rssi = randn()+signal_map[index[1],index[2]]
  #     signals[i] = RssiRecord(rssi,1,1.)
  #   end
  #   return path,signals
  # end

  # function rssi_distribution(rssi::Float64,signal_map::Array{Float64},spaceInfo::SpaceInfo)
  #   dist = scaled_normpdf(rssi,signal_map,.01,100.)
  #   return dist#/sum(dist)
  # end

  function rssiLogDist(rssi::Float64, # current rssi value
                      signalMap::Array{Float64} # unrolled coverage map
                      )
    # use coverage map values as mean, rssi as value and
    # rssiSigma as stdev for gaussian distribution
    dist = lognormpdf(rssi,signalMap,rssiSigma)
    # log probabilities are normalized
    # println(dist)
    # maximum = log(sum(exp.(dist)))
    # val = exp.(dist)
    return dist
    # return dist-maximum
  end

  # function transition_distribution(state,v,spaceInfo,dt = 1.,boosting = 1e100)
  #   location = grid2coord(fold_index(state,spaceInfo),spaceInfo.grid_size)
  #   dx = collect(0:1:(spaceInfo.grid[1]-1))*spaceInfo.grid_size+spaceInfo.grid_size/2
  #   dy = collect(0:1:(spaceInfo.grid[2]-1))*spaceInfo.grid_size+spaceInfo.grid_size/2
  #   x = scaled_normpdf(dx,location[1]+v[1]*dt,dt^2/2,50.)
  #   y = scaled_normpdf(dy,location[2]+v[2]*dt,dt^2/2,50.)'
  #   map = x*y
  #   ax = (location[1]+v[1]*dt-dx)*2/dt^2
  #   ay = (location[2]+v[2]*dt-dy)*2/dt^2
  #   vx = ones(Float64,spaceInfo.grid[2])*(v[1]+ax*dt)'
  #   vy = (ones(Float64,spaceInfo.grid[1])*(v[2]+ay*dt)')'
  #   # check velocity
  #   boosting_coefficient = boosting/sum(map)
  #   return map*boosting_coefficient,[vx[:] vy[:]]
  # end

  function transLogDist(state::Int, # previous state. defines mean value
                        v::Array{Float64}, # velocity defines mean value
                        plan,
                        dt::Float64)
    g_size = 1. # only one grid resolution is currently supported
    s_size = plan.sp2d_size
    limits = plan.limits
    # println("fold_index = $(fold_index(state,plan))")
    loc = grid2coord(fold_index(state,plan),plan)
    # create strided distance vector for further calculations
    # create these once outside this function
    dx = collect( limits[1,1]:g_size:(limits[1,2]-1) )  + g_size/2
    dy = collect( limits[2,1]:g_size:(limits[2,2]-1) )  + g_size/2
    # need support for rectangular space
    x = lognormpdf(dx,loc[1]+v[1]*dt,dt/1.414213562373095)
    y = lognormpdf(dy,loc[2]+v[2]*dt,dt/1.414213562373095)

    x1 = normalpdf(dx,loc[1]+v[1]*dt,dt/1.414213562373095)
    y1 = normalpdf(dy,loc[2]+v[2]*dt,dt/1.414213562373095)

    map1 = broadcast(*,x1,y1')
    log.(map1)


    # x = transProb(dx, loc[1], v[1], dt)
    # y = transProb(dy, loc[2], v[2], dt)
    # println("X = $(x)")
    # println("Y = $(y)")
    map = broadcast(+,x,y')
    # map = (x*ones(Float64,1,length(y)))+(ones(Float64,length(x),1)*y)
    ax = (dx-loc[1]-v[1]*dt)*2/dt^2
    ay = (dy-loc[2]-v[2]*dt)*2/dt^2
    #---------------------- println("The aceleration by x=$(ax) | y=$(ay)")
    vx = (ones(Float64,s_size[2])*(v[1]+ax*dt)')'
    vy = ones(Float64,s_size[1])*(v[2]+ay*dt)'
    #---------------------- println("The velocity by x=$(vx) | y=$(vy)")
    # check velocity
    sum(exp.(map))
    lpdfmax = maximum(map)            #deleted normolize?
    return map-lpdfmax,[vx[:] vy[:]]
  end

  function transProb(x, s, v, t, penalty = 1, acc_var = .30)

       # Calculate probability distribution of mobility model with velocity penalty
       # x : coordinate
       # v : initial velocity
       # t : time interval
       # penalty : number from 0 to inf, hom much penalty is assigned to the
       #             current value of velocity
       # acc_var : given that the acceleration has normal distribution, specifies
       #             acceleration variance

       v_abs = abs(v)
       # Parameters are based on the assumption that the pedastrian is most likely to
       # travel with the same speed as before
       mean = s + v * t
       # The value for variance is based on the motion equation
       hig_pen_var = acc_var * t^2 / (1 + (v>0) * v_abs * penalty) / 2
       low_pen_var = acc_var * t^2 / (1 + (v<0) * v_abs * penalty) / 2

       p_low_lb = mean - 5 * low_pen_var
       p_hig_ub = mean + 5 * hig_pen_var

       pu = TruncatedNormal(mean, hig_pen_var, p_low_lb, p_hig_ub)
       # pl = TruncatedNormal(mean, low_pen_var, p_low_lb, mean)

       # p_low = pdf.(pl, x);
       p_hig = pdf.(pu, x)
      #  println(sum(p_low)," ",sum(p_hig))

       # p_low_sum = sum(p_low)
       p_hig_sum = sum(p_hig)

       # if p_low_sum+p_hig_sum==0.  #x[end] < p_low_lb || x[1] > p_hig_ub
       if p_hig_sum==0.
          # if all probabilities are zero
          # println("\n All probabilities are zero")
          return p_hig-200.
       elseif  p_hig_sum > 0.
          # mean value is less than provided x
          p_tot = p_hig

       # elseif p_hig_sum == 0. && p_low_sum > 0.
       #    # mean value is greater than provided x
       #    p_tot = p_low
       # else
       #    nze = find(p_hig)[1] # find first non-zero
       #    factor = p_low[nze]/p_hig[nze]
       #    p_low[nze] = 0
       #    p_hig *= factor
       #    p_tot = (p_low+p_hig)
       end
       total_factor = 1/sum(p_tot)
       p_tot *= total_factor
       return map(p->if p>0 log10(p) else -200. end, p_tot)
    end


  # ==== Needs improvement ====
  function path_to_rssi(path,signal_map,plan, dt = 1.)
    #   this function needs to be adapted for the use with mapPlan
    g_size = 1.
    signalInfo = Array(RssiRecord,size(path,1))
    rssi = zeros(Float64,size(path,1),1)
    index_path = zeros(Float64,size(path,1),2)
    ipu = zeros(Float64,size(path,1),1)
    v1 = zeros(Float64,6,2)
    a1 = zeros(Float64,6,2)
    for i in 1:1:size(path,1)
      index = coord2grid(path[i,:],plan)
      rssi[i] = randn()*rssiSigma+signal_map[1][index[1],index[2]]
      index_path[i,:] = index
      ipu[i] = unfold_index(index,plan)

      signalInfo[i] = RssiRecord(rssi[i],1,1.)

      if i<size(path,1)
        a1[i,:] = (path[i+1,:]-path[i,:]-v1[i,:]*dt)*2/dt^2
        v1[i+1,:] = v1[i,:] + a1[i,:]*dt
      end

    end
    return signalInfo,index_path,ipu, v1
  end


  # function combine_prob(p_rssi,step,boosting = 1e100)
  #   prob = p_rssi.*step
  #   boosting_coefficient = boosting/sum(prob)
  #   return prob[:]*boosting_coefficient
  # end

  function combine_logprob(p_rssi::Array{Float64},step::Array{Float64})
    # println("Before combine_prob: p_rssi=$(p_rssi) | step=$(step)")
    prob = p_rssi+step
    # println("After combine_prob: prob=$(prob)")
    return prob[:] # should do normalization?
  end

  function combine_trellis(trellis_part,velocity_part,p_transition,v_transition,path_back,state)
    # vectorization possible?
    # ----------------------println("before combine_trellis: trellis_part=$(trellis_part) | velocity_part=$(velocity_part) | p_transition=$(p_transition) | v_transition=$(v_transition) | path_back=$(path_back) | state=$(state)")
    for ss in 1:1:length(trellis_part)
      if p_transition[ss] > trellis_part[ss]# maybe >=
        trellis_part[ss] = p_transition[ss]
        velocity_part[ss,:] = v_transition[ss,:]
        path_back[ss] = state
      elseif p_transition[ss] == trellis_part[ss]
        println("p_transition = $(p_transition[ss]) and Trellis = $(trellis_part[ss])")
      end
    end
    return trellis_part,velocity_part,path_back
  end

  # function combine_logtrellis(temp_trellis::Array{Float64},temp_velocity::Array{Float64})
  #   g_size = size(temp_trellis,1)
  #   trellis = zeros(Float64,g_size)
  #   velocity = zeros(Float64,g_size,2)
  #   path_back = zeros(Int32,g_size)
  #   (~,maxind) = findmax(temp_trellis,1)
  #
  #   trellis = temp_trellis[maxind]
  #   velocity = temp_velocity[maxind]
  #   temp = ceil(maind/g_size)
  #   path_back = maxind - (temp-1)*g_size
  #
  #   return trellis,velocity,path_back
  # end

  function beam_search(data, K) # data = trellis[:,1]
    trellis = sort(data, rev=true)
    state = findin(data,trellis[1:K])
    return state
    # sequences = [[list(), 1.0]]
    # 	# walk over each step in sequence
    # 	for row in data:
    # 		all_candidates = list()
    # 		# expand each current candidate
    # 		for i in range(len(sequences)):
    # 			seq, score = sequences[i]
    # 			for j in range(len(row)):
    # 				candidate = [seq + [j], score * -log(row[j])]
    # 				all_candidates.append(candidate)
    # 		# order all candidates by score
    # 		ordered = sorted(all_candidates, key=lambda tup:tup[1])
    # 		# select k best
    # 		sequences = ordered[:k]
    # 	return sequences
  end

  # function greedy_decoder(data) # input : trellis[:,1] it's vector
	# # index for largest probability each row
  # # version for julia
  #  # mxval, mxindx = findmax(data, ndims(data))
  #  # ind2sub(size(a), vec(mxindx))[2]
  #  return indmax(data)
	#  # return [argmax(s) for s in data] #python version
  # end

  function estimate_path_viterbi(signals,#::Array{RssiRecord},
                                  signal_map,#::Array{Array{Float64}},
                                  plan;
                                  seed = [-1. -1.],
                                  testing = false,
                                  ip = [],
                                  ipu = [])

    # This function estimates the trajectory of a network user based on
    # the history of observed signal strength. Input arguments are
    #
    # signals - sequence of observed rssi that contains rssi, assosiated AP id,
    #           and the time difference from the previous sample
    # signal_map - contains signal strength matrices that cover all the map for
    #               all available APs. Each matrix in signal_map is a signal
    #               strength map for an AP with id
    # spaceInfo - contains information about space configuration such as map
    #             grid size, the size of a cell on a grid, etc.
    # seed - initial location guessed by the viterbi algorithm. Used for testing
    #         purposes. May be excluded in the future.
    # testing - testing flag. Used to indicate whether additional output is
    #           needed
    # ip - ???
    # ipu - ???
    #
    # Probabilities in this function are subject to overflow and undeflow.
    # In order to cope with this problem probability scaling is used


    g_size = prod(plan.sp2d_size) # total number of cells in the grid
    trellis = zeros(Float64,g_size,2)-1e10 # stores previous and current state probability
    velocity = zeros(Float64,g_size,2,2) # stores previous and current state velocity in xy dimensions
    path_back = zeros(Int32,g_size,length(signals)) # stores feasible paths


    # initialize state probabilities based on seed location, or assume
    # equal probabilities
    init_step = zeros(Float64,tuple(plan.sp2d_size...))
    if seed==[-1. -1.]
      init_step = uniform_location(plan)
    else
      init_step = seed_location(seed,plan)#*1e100
    end


    # calculate log distribution. log allows to use sum instead of product
    # This is used to initialize probabilities in trellis
    c_rssi = signals[1].rssi
    c_ap_id = signals[1].ap
    # println("before rssiLogDist")
    p_rssi = rssiLogDist(c_rssi,
                          signal_map[c_ap_id]) # the value is boosted by exp(100)

    # Combines initial probabilities
    # println("combine_logprob")
    # @time
    trellis[:,1] = combine_logprob(p_rssi,init_step) # values sum to 1e100?


    # Start processing rssi signals
    # println("for i = 1:1:length(signals)-1")
    for i = 1:1:length(signals)-1
      # println("viterby Start processing rssi signals $(i)")
      # whenever there is no time difference between consecutive measurements,
      # assume that we did not change location
      dt = signals[i+1].t - signals[1].t
      if dt==0.
        for ii=1:length(path_back[:,i+1])
          path_back[ii,i+1] = ii
        end
        continue
      end

      # println("before rssiLogDist - Location probabilities based on rssi")
      # Location probabilities based on rssi
      c_rssi = signals[i+1].rssi
      c_ap_id = signals[i+1].ap
      p_rssi = rssiLogDist(c_rssi,
                            signal_map[c_ap_id])
      # if (i == 16 || i == 17)
      #   println("\n\n\nThe c_rssi of $(i)= $(c_rssi)")
      #   if (i == 16)
      #     println("\n\n\nThe p_rssi of $(i)= $(p_rssi[425])1749 and 1807 real = $(p_rssi[1807])")
      #   else
      #     println("\n\n\nThe p_rssi of $(i)= $(p_rssi[415])2298 and 1868 real =$(p_rssi[1868])")
      #   end
      # end
    #   if testing
    #     print("Maximum signal proability on the next step is at ",
    #           fold_index(indmax(p_rssi),plan),
    #           " and ",
    #           fold_index(ipu[i+1],plan),
    #           " is expected\n")
    #   end
      # println(trellis[:,1])
      # println("before for state in 1:1:g_size")

      # for state in 1:1:g_size
      for state in beam_search(trellis[:,1],50)
        # iterate over states in trellis for exhaustive optimum search
        # println("State = $(state)")
        # continue onlu if the log probability of the previous state is
        # more than -100. Remember scailing to 1e100 exists
        # beam_search(trellis[:,1],)
        # if (trellis[state,1]>-100)
          # Calculate transition probability based on time difference
          # and previous state

           # println("\n\n\n\n\n velocity= $(velocity[state,:,1])")
          p_transition,v_transition = transLogDist(state,
                                                  velocity[state,:,1],
                                                  plan,
                                                  dt)

        #   if testing
        #     if state == ipu[i]
        #       print("Maximum location proability on the next step is at ",
        #             fold_index(indmax(combine_logprob(p_rssi,p_transition)),plan),
        #             " and ",
        #             fold_index(ipu[i+1],spaceInfo),
        #             " is expected\n")
        #     end
        #   end
          # println("before combine_logprob")
          # if (i==16)
          #     if state == 425
          #     println("prob est_path for i p_rssi= $(p_rssi[state]) and p_transition= $(p_transition[state])")
          #   elseif state == 1807
          #     println("prob real_path for i p_rssi= $(p_rssi[state]) and p_transition= $(p_transition[state])")
          #   end
          # elseif i==17
          #     if state == 415
          #     println("prob est_path for i p_rssi= $(p_rssi[state]) and p_transition= $(p_transition[state])")
          #   elseif state == 1868
          #     println("prob real_path for i p_rssi= $(p_rssi[state]) and p_transition= $(p_transition[state])")
          #   end
          # end
          cp = combine_logprob(p_rssi,p_transition)
          # println("before combine_trellis")
          # if (i==16)
          #     if state == 425
          #       println("prob est_path for i cp= $(cp[state]) and trellis= $(trellis[state,1])")
          #   elseif state == 1807
          #     println("prob real_path for i cp= $(cp[state]) and trellis= $(trellis[state,1])")
          #   end
          # elseif i==17
          #     if state == 415
          #       println("prob est_path for i cp= $(cp[state]) and trellis= $(trellis[state,1])")
          #   elseif state == 1868
          #     println("prob real_path for i cp= $(cp[state]) and trellis= $(trellis[state,1])")
          #   end
          # end

          # slicing creates new arrays -> performance issue
          trellis[:,2],velocity[:,:,2],path_back[:,i+1] = combine_trellis(trellis[:,2],
                                                                          velocity[:,:,2],
                                                                          cp+trellis[state,1],
                                                                          v_transition,
                                                                          path_back[:,i+1],
                                                                          state)
          # if (i==16)
          #     if state == 425
          #       println("prob est_path for i trellis2= $(trellis[state,2])")
          #   elseif state == 1807
          #     println("prob real_path for i trellis2= $(trellis[state,2])")
          #   end
          # elseif i==17
          #     if state == 415
          #       println("prob est_path for i trellis2= $(trellis[state,2])")
          #   elseif state == 1868
          #     println("prob real_path for i  trellis2= $(trellis[state,2])")
          #   end
          # end

          # println(trellis[1,:])
          # println("The result after combine_trellis: $(path_back[:,i+1])")
          # ----------------------println("The result after combine_trellis: trellis2=$(trellis[:,2]) | velocity=$(velocity[:,:,2]) | $(path_back[:,i+1])")
          # println("------------------------------------")

        # end

      end
      # if (i == 16 || i == 17)
        # println("\n\n\nThe p_transition of $(i)= $(p_transition)")
        # println("\n\n\nThe v_transition of $(i)= $(v_transition)")
        # println("\n\n\nThe cp of $(i)= $(cp)")
        # println("\n\n\nThe v_transition of $(i)= $(v_transition)")

        # println("\n\n\nThe path_back[:,i] of $(i)= $(path_back[:,i])")
        # println("\n\n\nThe velocity[:,:,2] of $(i)= $(velocity[:,:,2])")
        # println("\n\n\nThe trellis[:,2] of $(i)= $(trellis[:,2])")

      # end
      # print("\n")
      # trellis[:,i+1],velocity[:,:,i+1],path_back[:,i+1] = combine_logtrellis(temp_trellis,temp_velocity)
      # boosting = 1e100
      # boosting_coefficient = boosting/sum(trellis[:,i+1])
      # trellis[:,i+1] = trellis[:,i+1]*boosting_coefficient

      # find the maximum probability value
      # ----------------------println("The second column or row trellis = $(trellis[:,2])")

      lpdfmax = maximum(trellis[:,2])
      # println("\n\n\nThe maximum value of trellis = $(lpdfmax) in o =$(i)")
      # println("\n\n\nindex of The maximum value trellis= $(find(a->a==lpdfmax, trellis[:,2]))")
      # # scale log probabilities? slicing creates new arrays -> performance issue
      trellis[:,2] = trellis[:,2]-lpdfmax
      # println("scale log probabilities trellis = $(trellis[:,2])")
      # current state probabilities and velocities assigned to the previous
      # step for the next iteration
      trellis[:,1] = trellis[:,2]
      velocity[:,:,1] = velocity[:,:,2]
      trellis[:,2] = zeros(Float64,g_size)-1e10 # creating new array -> performance issues
      velocity[:,:,2] = zeros(Float64,g_size,2) # creating new array -> performance issues
      # toc()
    end
    # println("Allocate memory for estimated path")
    # Allocate memory for estimated path
    estimated_path = zeros(Float64,length(signals),2)
    # Find most likely path
    step_back = indmax(trellis[:,1])
    # println("step_back= $(step_back)")
    # println("another way step_back= $(findmax(trellis[:,1]))")
    # A_max = maximum(trellis[:,1])
    # println("ind of all maxima= $(find(a->a==A_max, trellis[:,1]))")
    # print("Choosen end max is at $(step_back)\n")
    # estimated_path[end,:] = grid2coord(fold_index(location,spaceInfo),spaceInfo.grid_size)
    # step_back = path_back[location,end]
    steps = []
    # Propagate back to find the most likely path
    for step_ind in length(signals):-1:1
      # println("path-back= $(path_back) |step_ind= $(step_ind)")
      # println("step_back= $(step_back) |step_ind= $(step_ind)")
      push!(steps,step_back)
      estimated_path[step_ind,:] = grid2coord(
                                      fold_index(step_back,plan),
                                                  plan)
      step_back = path_back[step_back,step_ind]

    end


    # return trellis,velocity
    return estimated_path,trellis,steps
  end
end
