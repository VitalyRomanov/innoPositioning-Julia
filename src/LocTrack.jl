module LocTrack

  # dt = 1
  using MapPlan
  using Base.Test
  using Distributions
  # using Plots

  const rssiSigma = 5 # rssi noise stdev
  const sqrt_2 = sqrt(2)
  const sqrt_2pi = sqrt(2*pi)
  const eps_log = 1.e-200
  const log_eps_log = log(eps_log) # const

  export getSpaceInfo,estimate_path_viterbi,RssiRecord,pmodel_space,path_to_rssi,path_generation

  type RssiRecord
    rssi::Float64
    ap::Int
    t::Float64
  end

  type scaledPlan
      limits::Array{Int}
      gcell_size::Float64
      sp2d_size::Array{Int}
  end

  struct transLogDist_temp
      x_grid
      y_grid
      map_brc
      px
      py
      ax
      ay
      vx
      vy
      vx_brc
      vy_brc
      v_brc
      vx_ones
      vy_ones
  end

  transLogDist_temp(plan) = begin

      g_size = plan.gcell_size
      s_size = plan.sp2d_size
      limits = plan.limits

      x_grid = collect( limits[1,1]:g_size:(limits[1,2]-1) ) + g_size/2
      y_grid = (collect( limits[2,1]:g_size:(limits[2,2]-1) ) + g_size/2)'

      s_x = length(x_grid)
      s_y = length(y_grid)

      map_brc = zeros(s_x, s_y)
      vx_brc = zeros(s_x, s_y)
      vy_brc = zeros(s_x, s_y)

      vx_ones = zeros(1, s_y)
      vy_ones = zeros(s_x, 1)

      p_rssi = zeros(s_x, s_y)
      p_trans = zeros(s_x, s_y)
      combined = zeros(s_x, s_y)

      px = zeros(size(x_grid))
      py = zeros(size(y_grid))

      ax = zeros(size(x_grid))
      ay = zeros(size(y_grid))

      vx = zeros(size(x_grid))
      vy = zeros(size(y_grid))

      v_brc = zeros(s_x, s_y, 2)

      return transLogDist_temp(
          x_grid, y_grid,
          map_brc,
          px, py,
          ax, ay,
          vx, vy,
          vx_brc, vy_brc, v_brc,
          vx_ones, vy_ones
      )
  end

  struct viterbi_temp
      p_rssi
      combined
  end

  viterbi_temp(size) = begin
      p_rssi = zeros(size[1],size[2])
      combined = zeros(size[1],size[2])
      return viterbi_temp(p_rssi, combined)
  end

  struct transProbCoord_temp
      mean
      p_low_lb
      p_hig_ub
      p_lower
      p_upper
  end

  struct TestPlan
      limits
      gcell_size
      sp2d_size
  end

  ##### Auxiliary functions


  # log of Normal PDF
  lognormpdf(x,m,s) = -((x.-m)./s).^2/2-log(sqrt_2pi.*s)

  # coord2grid takes location on the map and converts it to grid coords
  coord2grid(coord,plan) = Int.(ceil.( (coord-plan.limits[1:2,1]) / plan.gcell_size ))

  # grid2coord takes grid coords and converts it to location on the map
  grid2coord(grid_loc,plan) = grid_loc*plan.gcell_size + plan.limits[1:2,1] - plan.gcell_size/2 # .5 is the half of sector size of 1 m2

  # unfold_index taked 2d grid location and does mapping from 2d to 1d index
  unfold_index(grid_loc,plan) = grid_loc[1]+(grid_loc[2]-1)*plan.sp2d_size[1]

  # fold_index taked 1d index and does mapping from 1d to 2d grid location
  fold_index(index,plan) = begin
      ii = [0,0]
      ii[2] = ceil(index/plan.sp2d_size[1])
      ii[1] = index[1] - (ii[2]-1)* plan.sp2d_size[1]
      ii
  end

  function seed_location(seed,plan)             # +1 sec in running time viterbi
    # creates distrubution over space with all probability mass in the position
    # specified by seed
      index = coord2grid(seed,plan)
      seeded_location = zeros(Float64,tuple(plan.sp2d_size...))+eps_log
      seeded_location[index[1],index[2]] = 1.
      seeded_location .=log.(seeded_location)
      return seeded_location
  end

  uniform_location(plan) = begin
    #   creates a distribution over space with all positions beimg equally likely
      spacesize = prod(plan.sp2d_size)
      ones(Float64,tuple(plan.sp2d_size...))/spacesize
  end

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

  function path_generation(distribution,
                            path_lenght,
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
   r_sigma = 1 #instead 0.1
   t = dt
   path = zeros(Float64,path_lenght,2) # stores generated path here
   rssi = zeros(Float64,path_lenght,1)
   signals  = Array(RssiRecord,path_lenght)
   v = zeros(Float64,2,1)
   path[1,:] = seed

   index = coord2grid(path[1,:],plan)
   # print(index)
   # rssi[1] = randn()+signal_map[index[1],index[2]]
   rssi = signal_map[1][index[1],index[2]]
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
     rssi = signal_map[1][index[1],index[2]]
     signals[i] = RssiRecord(rssi,1,t)
   end
   return path,signals
  end



  # function rssiLogDist!(p_rssi,
  #                     rssi::Float64, # current rssi value
  #                     signalMap::Array{Float64},
  #                     sigma = 1.)
  #   # use coverage map values as mean, rssi as value and rssiSigma as stdev for gaussian distribution
  #   @. p_rssi = lognormpdf(rssi,signalMap,sigma)
  #   max_val = maximum(p_rssi)
  #   if isinf(max_val)
  #       p_rssi = 0.
  #   else
  #       p_rssi.-=max_val
  #   end
  #   # return p_rssi
  # end

  function rssiLogDist(rssi::Float64, # current rssi value
                      signalMap::Array{Float64},
                      sigma = 1.)
    # use coverage map values as mean, rssi as value and rssiSigma as stdev for gaussian distribution
    p_rssi = lognormpdf(rssi,signalMap,sigma)
    max_val = maximum(p_rssi)
    if isinf(max_val)
        p_rssi = 0.
    else
        p_rssi.-=max_val
    end
    # return p_rssi
  end

  function spaceInit(plan)

      tld_temp = transLogDist_temp(plan)
      vit_temp = viterbi_temp(plan.sp2d_size)

      return tld_temp, vit_temp
  end

  function transProbCoord(x, s, v, t, penalty = 10, acc_var =1)# 0.5487618458)
         v_abs = abs(v)
         mean = s + v * t
         hig_pen_var = acc_var * t^2 / (1 + (v>0) * v_abs * penalty) / 2
         low_pen_var = acc_var * t^2 / (1 + (v<0) * v_abs * penalty) / 2
         p_low_lb = mean - 5 * low_pen_var              #жестко ограничен никуда не двигается на рабочих данных
         p_hig_ub = mean + 5 * hig_pen_var
         pu = TruncatedNormal(mean, sqrt(hig_pen_var), p_low_lb, p_hig_ub)
         # pl = TruncatedNormal(mean, sqrt(low_pen_var), p_low_lb, mean)
         # pu = TruncatedNormal(mean, sqrt(hig_pen_var), mean, p_hig_ub)
         p_hig = pdf.(pu, x)
         p_hig_sum = sum(p_hig)
         if p_hig_sum==0.
             # println("p_hig_sum = 0 -> p_hig-200.")
            return p_hig + log_eps_log
         elseif  p_hig_sum > 0.
             # println("else p_tot")
            p_tot = p_hig
         end
         p_v = lognormpdf(v,0,2.1)
         p_v = exp(p_v)
         p_tot *= p_v
         return map(p->if p>0 log(p) else log(eps_log) end, p_tot)
  end


  # function transLogDist(state::Int, # previous state. defines mean value
  #                       v::Array{Float64}, # velocity defines mean value
  #                       plan,
  #                       dt::Float64,
  #                       dx,
  #                       dy,
  #                       s_size)
  #   g_size = plan.gcell_size # grid resolution
  #   s_size = plan.sp2d_size
  #   limits = plan.limits
  #
  #   loc = grid2coord(fold_index(state,plan),plan)
  #
  #   # create strided distance vector for further calculations
  #   # create these once outside this function
  #
  #   # need support for rectangular space
  #
  #   x = transProbCoord(dx,loc[1],v[1],dt)
  #   y = transProbCoord(dy,loc[2],v[2],dt)
  #   brc = broadcast(+,x,y')
  #   brc .= exp.(brc)
  #   brc ./= sum(brc)
  #   brc .= log.(brc)
  #
  #   ax = (dx - loc[1] - v[1] * dt) * 2 / dt^2
  #   ay = (dy - loc[2] - v[2] * dt) * 2 / dt^2
  #
  #   vx = broadcast(+,zeros(1,s_size[2]), v[1]+ax*dt)
  #   vy = broadcast(+,zeros(s_size[1],1), v[2]+ay'*dt)
  #
  #   return brc,[vx[:] vy[:]]
  # end

  function transLogDistFast(state::Int, # previous state. defines mean value
                        v::Array{Float64}, # velocity defines mean value
                        plan,
                        dt::Float64,
                        vars)

    loc = grid2coord(fold_index(state,plan),plan)

    dx = vars.x_grid
    dy = vars.y_grid

    brc = vars.map_brc
    vx_brc = @view vars.v_brc[:,:,1]
    vy_brc = @view vars.v_brc[:,:,2]

    px = vars.px
    py = vars.py

    ax = vars.ax
    ay = vars.ay

    vx = vars.vx
    vy = vars.vy

    @. px = lognormpdf(dx, loc[1] + v[1] * dt, dt / sqrt_2) + lognormpdf(v[1], 0, 1.5)
    @. py = lognormpdf(dy, loc[2] + v[2] * dt, dt / sqrt_2) + lognormpdf(v[1], 0, 1.5)

    # px = transProbCoord(dx,loc[1],v[1],dt)
    # py = transProbCoord(dy,loc[2],v[2],dt)

    broadcast!(+,brc,px,py)
    brc .= exp.(brc)
    brc ./= sum(brc)
    brc .= log.(brc)

    @. ax = (dx - loc[1] - v[1] * dt) * 2 / dt^2
    @. ay = (dy - loc[2] - v[2] * dt) * 2 / dt^2

    @. vx = v[1]+ax*dt
    @. vy = v[2]+ay*dt

    # println(size(vars.vx_brc), size(vars.vx_ones), size(v[1]+ax*dt))
    broadcast!(+, vx_brc, vars.vx_ones, vx)
    broadcast!(+, vy_brc, vars.vy_ones, vy)

    return brc, vars.v_brc
  end

  function combine_logprob!(prob, p_rssi, step)
    @. prob = p_rssi+step
    prob.-= maximum(prob)
    prob .= exp.(prob)
    prob./=sum(prob)
    prob .= log.(prob)
    # return prob[:]
  end

  function combine_trellis(trellis_part,velocity_part,p_transition,v_transition,path_back,state)
    # vectorization possible?
    for ss in 1:1:length(trellis_part)
      if p_transition[ss] > trellis_part[ss]
        trellis_part[ss] = p_transition[ss]
        velocity_part[ss,:] = v_transition[ss,:]
        path_back[ss] = state
      end
    end
    return trellis_part,velocity_part,path_back
  end

  function combine_trellis_fast!(trellis,velocity,p_transition,v_transition,path_back,t_next,state, plan)

    for ss = 1:length(p_transition)
      if p_transition[ss] > trellis[ss,2]
        trellis[ss,2] = p_transition[ss]
        ind = fold_index(ss,plan)
        velocity[ss,:,2] = @view v_transition[ind[1],ind[2],:]
        path_back[ss,t_next] = state
      end
    end
  end

  function weight_with_state!(cp, p_state)
      cp .+= p_state
  end

  function scaleSpace(plan,ssms,factor)
      newssms = Array{Array{Float64}}(0)
      new_size = Int.(ceil.( (plan.limits[1:2,2]-plan.limits[1:2,1]) / factor))
      for ssm in ssms
          newssm = zeros((new_size...))
          for x=1:size(newssm,1), y=1:size(newssm,2)
              x_top = min(size(ssm,1),x*factor-1)
              y_top = min(size(ssm,2),y*factor-1)
              newssm[x,y] = maximum(ssm[x*factor-factor+1:x_top,y*factor-factor+1:y_top])
          end
          push!(newssms,newssm)
      end

      newplan = scaledPlan(plan.limits,factor,new_size)
      return newplan,newssms
  end



  function beam_search(data, K) # data = trellis[:,1]
      states = sortperm(data, rev=true)
    # trellis = sort(data, rev=true)
    # state = findin(data,trellis[1:K])
    return states[1:K]
  end


  function estimate_path_viterbi(signals,
                                  o_signal_map,
                                  o_plan;
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


    plan,signal_map = scaleSpace(o_plan,o_signal_map,5)
    tdl_t, vtrb_t = spaceInit(plan)

    p_rssi = vtrb_t.p_rssi
    cp = vtrb_t.combined


    g_size = prod(plan.sp2d_size) # total number of cells in the grid
    trellis = zeros(Float64,g_size,2)-Inf # stores previous and current state probability
    velocity = zeros(Float64,g_size,2,2) # stores previous and current state velocity in xy dimensions
    path_back = zeros(Int32,g_size,length(signals)) # stores feasible paths
    s_size = plan.sp2d_size



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
    # c_rssi = signals[1].rssi
    # c_ap_id = signals[1].ap
    # rssiLogDist!(p_rssi, c_rssi, signal_map[c_ap_id])

    p_rssi[:] = 1.
    for sig in signals[1]
        p_rssi += rssiLogDist(sig.rssi, signal_map[sig.ap],rssiSigma)
    end

    # Combines initial probabilities
    combine_logprob!(cp, p_rssi, init_step)
    trellis[:,1] = cp[:]


    # Start processing rssi signals
    for i = 1:1:length(signals)-1
      # whenever there is no time difference between consecutive measurements,
      # assume that we did not change location
      dt = signals[i+1][1].t - signals[i][1].t
      # if dt==0.
      #   for ii=1:length(path_back[:,i+1])
      #     path_back[ii,i+1] = ii
      #   end
      #   continue
      # end

      # Location probabilities based on rssi
      # c_rssi = signals[i+1].rssi
      # c_ap_id = signals[i+1].ap
      # rssiLogDist!(p_rssi, c_rssi, signal_map[c_ap_id])

      # println(signals[i+1])
      p_rssi[:] = 1.
      for sig in signals[i+1]
          p_rssi += rssiLogDist(sig.rssi, signal_map[sig.ap],rssiSigma)
      end

        states = beam_search(trellis[:,1], 30)#div(length(trellis[:,1]),60))
        # println("t=",i, " ", length(states), " ", size(trellis[:,1]))

        for state in states
          # Calculate transition probability based on time difference
          # and previous state

          p_transition, v_transition = transLogDistFast(state,
                                                  velocity[state,:,1],
                                                  plan,
                                                  dt,
                                                  tdl_t)

          combine_logprob!(cp, p_rssi, p_transition)

          weight_with_state!(cp,trellis[state,1])
          combine_trellis_fast!(trellis,
                                velocity,
                                cp,
                                v_transition,
                                path_back,
                                i+1,
                                state,
                                plan)
      end
      trellis[:,1] = @view trellis[:,2]
      velocity[:,:,1] = @view velocity[:,:,2]
      trellis[:,2] = -Inf
      velocity[:,:,2] = 0.

      # ts = sort(trellis[:,1], rev=true)
      # plot(ts[1:20])
      # savefig("tdist/$(i+1).png")
    end
    # Allocate memory for estimated path
    estimated_path = zeros(Float64,length(signals),3)
    # Find most likely path
    step_back = indmax(trellis[:,1])

    # Propagate back to find the most likely path
    for step_ind in length(signals):-1:1
      estimated_path[step_ind,1:2] = grid2coord(fold_index(step_back,plan),plan)
      estimated_path[step_ind,3] = signals[step_ind][1].t
      step_back = path_back[step_back,step_ind]
    end

    return estimated_path
  end
end
