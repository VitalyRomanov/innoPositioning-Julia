module LocTrack

  # dt = 1

  const rssiSigma = 0.1 # rssi noise stdev

  export getSpaceInfo,estimate_path_viterbi,RssiRecord,PModelParam,pmodel_space,path_to_rssi,path_generation

  type RssiRecord
    rssi::Float64
    ap::Int
    dt::Float64
  end

  type SpaceInfo
    space::Array{Int}
    grid::Array{Int}
    grid_size::Float64
  end

  type PModelParam
    alpha::Float32
    d0::Float32
    p_d0::Float32
  end

  ##### Auxiliary functions

  # modify in accordance with performance guidelines
  # function normpdf(x,m,s)
  #   return exp(-((x-m)./s).^2/2)./(sqrt(2*pi)*s)
  # end
  #
  # function scaled_normpdf(x::Float64,m::Array{Float64},s::Float64,scale::Float64)
  #   return exp(scale-((x-m)./s).^2/2)./(2.5066282746310002*s) # sqrt(2*pi)
  # end

  # sqrt(2pi) = 2.5066282746310002
  lognormpdf(x,m,s)  = @. -((x-m)/s)^2/2-log(2.5066282746310002*s)

  # function lognormpdf(x::Float64,m::Array{Float64},s::Float64)
  #   return -((x-m)./s).^2/2-log(2.5066282746310002*s) # sqrt(2*pi)
  # end

  # function lognormpdf(x::Array{Float64},m::Float64,s::Float64)
  #   return -((x-m)./s).^2/2-log(2.5066282746310002*s) # sqrt(2*pi)
  # end

  function getSpaceInfo(space,grid_size=1.)
    grid = convert(Array{Int},floor(space/grid_size))
    return SpaceInfo(space,grid,grid_size)
  end

  function pmodel_space(ap_loc,spaceInfo::SpaceInfo,pmodel::PModelParam)
    # creates coverage map
    x = ones(Float64,spaceInfo.grid[1],1) * collect(0:1:(spaceInfo.grid[2]-1))' * spaceInfo.grid_size
    y = (ones(Float64,spaceInfo.grid[2],1) * collect(0:1:(spaceInfo.grid[1]-1))')' * spaceInfo.grid_size
    d = sqrt((x-ap_loc[1]).^2+(y-ap_loc[2]).^2)
    d[1,1] = .1
    value = pmodel.p_d0 - 10 * pmodel.alpha * log10(d/pmodel.d0)
    return value
  end

  function coordinates_to_grid(coord::Array{Float64},grid_size::Float64)
    index = [0 0]
    index[1] = convert(Int32,ceil(coord[1]/grid_size))
    index[2] = convert(Int32,ceil(coord[2]/grid_size))
    return index
  end

  function grid_to_coordinates(index::Array{Int},grid_size::Float64)
    return index*grid_size-grid_size/2
  end


  function unfold_index(index::Array{Int},spaceInfo::SpaceInfo)
    return index[1]+(index[2]-1)*(spaceInfo.grid[1])
  end

  function fold_index(index,spaceInfo)
    ii = zeros(Int64,2,1)
    ii[2] = ceil(index/spaceInfo.grid[1])
    ii[1] = index[1] - (ii[2]-1)* spaceInfo.grid[1]
    return ii
  end

  function seed_location(seed,spaceInfo)
    index = coordinates_to_grid(seed,spaceInfo.grid_size)
    seeded_location = zeros(Float64,spaceInfo.grid[1],spaceInfo.grid[2])
    seeded_location[index[1],index[2]] = 1.
    return seeded_location
  end

  function uniform_location(spaceInfo)
    spacesize = spaceInfo.grid[1]*spaceInfo.grid[2]
    return ones(Float64,spaceInfo.grid[1],spaceInfo.grid[2])/spacesize
  end

  ####### Project functions

  function path_generation(seed,path_lenght,signal_map,spaceInfo, dt = 1.)
    # seed : initial location
    # path_lenght : --||--
    # signal_map : signal strength matrix
    # spaceInfo :
    # dt :
    path = zeros(Float64,path_lenght,2) # stores generated path here
    # rssi = zeros(Float64,path_lenght,1)
    signals  = Array(RssiRecord,path_lenght)
    v = zeros(Float64,2,1)
    path[1,:] = seed
    index = coordinates_to_grid(seed,spaceInfo.grid_size)
    # print(index)
    # rssi[1] = randn()+signal_map[index[1],index[2]]
    rssi = randn()+signal_map[index[1],index[2]]
    signals[1] = RssiRecord(rssi,1,1.)

    for i in 2:1:path_lenght
      valid_sample = false
      s = seed
      a = seed
      while !valid_sample
        a = randn(2,1)
        s = path[i-1,:] + v*dt + a*dt^2/2
        if (s[1]>0) && (s[2]>0) && (s[1]<spaceInfo.space[1,1]) && (s[2]<spaceInfo.space[2,1])
          valid_sample = true
        end
      end
      path[i,:] = s
      v += a*dt
      index = coordinates_to_grid(path[i,:],spaceInfo.grid_size)
      rssi = randn()+signal_map[index[1],index[2]]
      signals[i] = RssiRecord(rssi,1,1.)
    end
    return path,signals
  end

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
    return dist-maximum(dist)
  end

  # function transition_distribution(state,v,spaceInfo,dt = 1.,boosting = 1e100)
  #   location = grid_to_coordinates(fold_index(state,spaceInfo),spaceInfo.grid_size)
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
                        si::SpaceInfo,
                        dt::Float64)
    location = grid_to_coordinates(fold_index(state,si),si.grid_size) # location in x,y defined by state
    dx = collect(0:1:(si.grid[1]-1))*si.grid_size+si.grid_size/2
    dy = collect(0:1:(si.grid[2]-1))*si.grid_size+si.grid_size/2
    # need support for rectangular space
    x = lognormpdf(dx,location[1]+v[1]*dt,dt^2/2)
    y = lognormpdf(dy,location[2]+v[2]*dt,dt^2/2)'
    map = (x*ones(Float64,1,length(y)))+(ones(Float64,length(x),1)*y)
    ax = (dx-location[1]-v[1]*dt)*2/dt^2
    ay = (dy-location[2]-v[2]*dt)*2/dt^2
    vx = (ones(Float64,si.grid[2])*(v[1]+ax*dt)')'
    vy = ones(Float64,si.grid[1])*(v[2]+ay*dt)'
    # check velocity
    lpdfmax = maximum(map)
    return map-lpdfmax,[vx[:] vy[:]]
  end

  function path_to_rssi(path,signal_map,spaceInfo, dt = 1.)
    signalInfo = Array(RssiRecord,size(path,1))
    rssi = zeros(Float64,size(path,1),1)
    index_path = zeros(Float64,size(path,1),2)
    ipu = zeros(Float64,size(path,1),1)
    v1 = zeros(Float64,6,2)
    a1 = zeros(Float64,6,2)
    for i in 1:1:size(path,1)
      index = coordinates_to_grid(path[i,:],spaceInfo.grid_size)
      rssi[i] = randn()*rssiSigma+signal_map[index[1],index[2]]
      index_path[i,:] = index
      ipu[i] = unfold_index(index,spaceInfo)

      signalInfo[i] = RssiRecord(rssi[i],1,1.)

      if i<size(path,1)
        a1[i,:] = (path[i+1,:]-path[i,:]-v1[i,:]*dt)*2/dt^2
        v1[i+1,:] = v1[i,:] + a1[i,:]*dt
      end

    end
    return signalInfo,index_path,ipu, v1
  end


  function combine_prob(p_rssi,step,boosting = 1e100)
    prob = p_rssi.*step
    boosting_coefficient = boosting/sum(prob)
    return prob[:]*boosting_coefficient
  end

  function combine_logprob(p_rssi::Array{Float64},step::Array{Float64})
    prob = p_rssi+step
    return prob[:]
  end

  function combine_trellis(trellis_part,velocity_part,p_transition,v_transition,path_back,state)
    for ss in 1:1:length(trellis_part)
      if p_transition[ss] > trellis_part[ss]
        trellis_part[ss] = p_transition[ss]
        velocity_part[ss,:] = v_transition[ss,:]
        path_back[ss] = state
      end
    end
    return trellis_part,velocity_part,path_back
  end

  function combine_logtrellis(temp_trellis::Array{Float64},temp_velocity::Array{Float64})
    grid_space_size = size(temp_trellis,1)
    trellis = zeros(Float64,grid_space_size)
    velocity = zeros(Float64,grid_space_size,2)
    path_back = zeros(Int32,grid_space_size)
    (~,maxind) = findmax(temp_trellis,1)

    trellis = temp_trellis[maxind]
    velocity = temp_velocity[maxind]
    temp = ceil(maind/grid_space_size)
    path_back = maxind - (temp-1)*grid_space_size

    return trellis,velocity,path_back
  end




  function estimate_path_viterbi(signals,
                                  signal_map,
                                  spaceInfo;
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


    grid_space_size = spaceInfo.grid[1]*spaceInfo.grid[2] # total number of cells in the grid
    # trellis = zeros(Float64,grid_space_size,length(signals))-1e10
    # velocity = zeros(Float64,grid_space_size,2,length(signals))
    trellis = zeros(Float64,grid_space_size,2)-1e10 # stores previous and current state probability
    velocity = zeros(Float64,grid_space_size,2,2) # stores previous and current state velocity in xy dimensions
    path_back = zeros(Int32,grid_space_size,length(signals)) # stores feasible paths


    # initialize state probabilities based on seed location, or assume
    # equal probabilities
    init_step = zeros(Float64,spaceInfo.grid[1],spaceInfo.grid[2])
    if seed==[-1. -1.] # first step
      init_step = uniform_location(spaceInfo) #set probability 1/N*M for each cell
    else
      init_step = seed_location(seed,spaceInfo)#*1e100
      if testing
        print("Initial location seeded at, ",
              fold_index(indmax(init_step),spaceInfo),
              " and ",
              fold_index(ipu[1],spaceInfo),
              " is expected\n")
      end
    end


    # calculate log distribution. log allows to use sum instead of product
    # This is used to initialize probabilities in trellis
    p_rssi = rssiLogDist(signals[1].rssi,
                          signal_map[signals[1].ap]) # the value is normolized

    # if testing
    #   print("Maximum signal proability is at ",
    #           fold_index(indmax(p_rssi),spaceInfo),
    #           " and ",
    #           fold_index(ipu[1],spaceInfo),
    #           " is expected \n")
    # end

    # Combines initial probabilities
    trellis[:,1] = combine_logprob(p_rssi,init_step) # values sum to 1e100


    # Start processing rssi signals
    for i = 1:1:length(signals)-1

      # whenever there is no time difference between consecutive measurements,
      # assume that we did not change location
      if signals[i+1].dt==0.
        for ii=1:length(path_back[:,i+1])
          path_back[ii,i+1] = ii
        end
        continue
      end

      # tic()
      # println("Step $i, dt = $(signals[i+1].dt)")

      # Location probabilities based on rssi
      p_rssi = rssiLogDist(signals[i+1].rssi,
                            signal_map[signals[i+1].ap])

      if testing
        print("Maximum signal proability on the next step is at ",
              fold_index(indmax(p_rssi),spaceInfo),
              " and ",
              fold_index(ipu[i+1],spaceInfo),
              " is expected\n")
      end


      for state in 1:1:grid_space_size
        # iterate over states in trellis for exhaustive optimum search

        # continue onlu if the log probability of the previous state is
        # more than -100. Remember scailing to 1e100 exists
        if (trellis[state,1]>-100)
          # Calculate transition probability based on time difference
          # and previous state
          p_transition,v_transition = transLogDist(state,
                                                  velocity[state,:,1],
                                                  spaceInfo,
                                                  signals[i+1].dt)

          if testing
            if state == ipu[i]
              print("Maximum location proability on the next step is at ",
                    fold_index(indmax(combine_logprob(p_rssi,p_transition)),spaceInfo),
                    " and ",
                    fold_index(ipu[i+1],spaceInfo),
                    " is expected\n")
            end
          end

          cp = combine_logprob(p_rssi,p_transition)

          # slicing creates new arrays -> performance issue
          trellis[:,2],velocity[:,:,2],path_back[:,i+1] = combine_trellis(trellis[:,2],
                                                                          velocity[:,:,2],
                                                                          cp+trellis[state,1],
                                                                          v_transition,
                                                                          path_back[:,i+1],
                                                                          state)
        end

      end
      print("\n")
      # trellis[:,i+1],velocity[:,:,i+1],path_back[:,i+1] = combine_logtrellis(temp_trellis,temp_velocity)
      # boosting = 1e100
      # boosting_coefficient = boosting/sum(trellis[:,i+1])
      # trellis[:,i+1] = trellis[:,i+1]*boosting_coefficient

      # find the maximum probability value
      lpdfmax = maximum(trellis[:,2])
      # scale log probabilities? slicing creates new arrays -> performance issue
      trellis[:,2] = trellis[:,2]-lpdfmax
      # current state probabilities and velocities assigned to the previous
      # step for the next iteration
      trellis[:,1] = trellis[:,2]
      velocity[:,:,1] = velocity[:,:,2]
      trellis[:,2] = zeros(Float64,grid_space_size)-1e10 # creating new array -> performance issues
      velocity[:,:,2] = zeros(Float64,grid_space_size,2) # creating new array -> performance issues
      # toc()
    end

    # Allocate memory for estimated path
    estimated_path = zeros(Float64,length(signals),2)
    # Find most likely path
    step_back = indmax(trellis[:,1])
    # print("Choosen end max is at $(step_back)\n")
    # estimated_path[end,:] = grid_to_coordinates(fold_index(location,spaceInfo),spaceInfo.grid_size)
    # step_back = path_back[location,end]

    # Propagate back to find the most likely path
    for step_ind in length(signals):-1:1
      estimated_path[step_ind,:] = grid_to_coordinates(
                                      fold_index(step_back,spaceInfo),
                                                  spaceInfo.grid_size
                                      )+spaceInfo.space[:,1]
      step_back = path_back[step_back,step_ind]
    end

    # return trellis,velocity
    return estimated_path,trellis
  end
end
