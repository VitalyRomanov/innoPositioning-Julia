module ClientDataProc
  using PyPlot
  using LocTrack
  using DataFrames

  export ip_to_labels!,filter_rssi_values!,change_rssi_scale!,remove_milliseconds,get_rssi_records,time_to_ms,filter_and_process

  function filter_and_process(input_file,id)
    client = DataFrames.readtable(input_file)
    client = ip_to_labels!(client)
    client = filter_rssi_values!(client,[0,1])
    client = change_rssi_scale!(client)
    client = remove_milliseconds(client)
    client[:,5] = deepcopy(client[:,4])
    client = ClientDataProc.time_to_ms(client)
    writetable("/Users/LTV/Documents/rssi_data/client_data/client_$(id).csv",client)
  end


  function get_rssi_records(client_data,number_of_records = -1)
    data_size = size(client_data,1)

    if number_of_records == -1
      stop_point = 1
    else
      stop_point = data_size-number_of_records
    end

    rssi_client_data = Array(RssiRecord,data_size)
    for i in data_size:-1:stop_point
      rssi_client_data[data_size-i+1] = RssiRecord(Float64(client_data[i,2]),client_data[i,3],client_data[i,4]/1000.)
    end

    timestamps = Array(client_data[:,5])

    if number_of_records == -1
      return rssi_client_data,timestamps
    else
      return rssi_client_data[1:number_of_records],timestamps[1:number_of_records]
    end

  end

  function get_ip_labels()
    ip_labeling = Dict()
    ips = ["11.0.0.1","11.0.0.2","11.0.0.117","11.0.0.21","11.0.0.22","11.0.0.25","11.0.0.26","11.0.0.23","11.0.0.24","11.0.0.10","11.0.0.11","11.0.0.14","11.0.0.13","11.0.0.18","11.0.0.19","11.0.0.12","11.0.0.17","11.0.0.15","11.0.0.7","10.200.2.1","10.200.0.1","10.200.1.1","10.200.2.2","10.200.0.2","10.200.0.177","10.200.2.3","10.200.1.3","10.200.0.21","10.200.0.22","10.200.2.5","10.200.0.25","10.200.0.26","10.200.2.7","10.200.0.23","10.200.0.24","10.200.1.5","10.200.2.8","10.200.0.10","10.200.0.11","10.200.2.9","10.200.0.14","10.200.0.13","10.200.1.6","10.200.2.15","10.200.0.18","10.200.0.19","10.200.2.119","10.200.0.12","10.200.0.17","10.200.0.15","10.200.0.7"]
    ids = [1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,9,10,1,1,1,2,2,2,3,3,3,3,4,4,4,5,5,5,5,6,6,6,7,7,7,7,8,8,8,9,9,9,9,10]
    # ip_labeling.keys = ips
    for i in 1:length(ips)
      ip_labeling[ips[i]] = ids[i]
    end
    return ip_labeling
  end

  function ip_to_labels!(client_data)
    i = 1
    client_record_size = size(client_data,1)
    ip_map = get_ip_labels()
    while i<=client_record_size
      ap_label = get(ip_map, client_data[i,:rssilog_antenna], -1)
      if ap_label == -1
        client_data = append!(client_data[1:i-1,:],client_data[i+1:end,:])
        i-=1
      else
        client_data[i,:rssilog_antenna] = "$(ap_label)"
      end
      i+=1
      client_record_size = size(client_data,1)
    end
    return client_data
  end

  function filter_rssi_values!(client_data,values)
    i = 1
    client_record_size = size(client_data,1)
    while i<=client_record_size
      if client_data[i,2] in values
        client_data = append!(client_data[1:i-1,:],client_data[i+1:end,:])
        i-=1
      end
      i+=1
      client_record_size = size(client_data,1)
    end
    return client_data
  end

  function select_ap(client_data,ap)
    i = 1
    client_record_size = size(client_data,1)
    while i<=client_record_size
      if client_data[i,3] != ap
        client_data = append!(client_data[1:i-1,:],client_data[i+1:end,:])
        i-=1
      end
      i+=1
      client_record_size = size(client_data,1)
    end
    return client_data
  end

  function plot_hist_by_ap(client_data,column)
    # select = deepcopy(client_data)
    for i in 1:10
      select = deepcopy(client_data)
      select = select_ap(select,"$i")
      h = PyPlot.plt[:hist](select[:,column],100)
      savefig("$(i)_ap.png")
      close()
    end
  end


  function change_rssi_scale!(client_data)
    for i in 1:size(client_data,1)
      client_data[i,2] = round(client_data[i,2]/2)
    end
    return client_data
  end

  function remove_milliseconds(client_data)
    ts = client_data[:,:rssilog_timestamp2]

    get_date(time) = split(time,".")[1]
    date2unix(datetime) = Dates.datetime2unix(Dates.DateTime(datetime,"yyyy-mm-dd HH:MM:SS"))
    convert(x) = begin
        try
            date2unix(get_date(x))
        catch
            println(x)
        end
    end

    ts = map(convert, ts)
    # ts = map(x->Dates.datetime2unix(Dates.DateTime(split(x,".")[1],"yyyy-mm-dd HH:MM:SS")), ts)

    client_data[:,:rssilog_timestamp2] = ts

    # client_data[i,:rssilog_timestamp2] = ts
    #
    # for i in 1:size(client_data,1)
    #   timestamp = split(client_data[i,:rssilog_timestamp2],".")
    #   client_data[i,:rssilog_timestamp2] = Dates.datetime2unix(Dates.DateTime(timestamp[1],"yyyy-mm-dd HH:MM:SS"))
    # end
    return client_data
  end

  function time_to_ms(client_data)
    timestamps = Array(DateTime,0)
    for i in 1:size(client_data,1)
      push!(timestamps,Dates.DateTime(client_data[i,4]))
    end
    tm = minimum(timestamps)
    timestamps = timestamps - tm
    for i in 1:size(client_data,1)-1
      # client_data[i,4] = split("$(timestamps[i])"," ")[1]
      client_data[i,4] = split("$(timestamps[i]-timestamps[i+1])"," ")[1]
    end
    client_data[size(client_data,1),4] = "0"
    return client_data
  end


  function get_client_pos(data)
      ids = :session_sso
      clients = Dict()

      n_records = size(data,1)
      for i = 1:n_records
          if i%1000 == 0
              print("\r$(i)/$(n_records) collecting positions")
          end
          client_id = data[i,ids]
          if client_id in keys(clients)
              push!(clients[client_id], i)
          else
              clients[client_id] = [i]
          end
      end

      return clients
  end


  function sift_clients(data)
      clients = get_client_pos(data)

      client_tables = Dict()
      println("")
      for (cl_i,cl) in enumerate(keys(clients))
          if cl_i % 100 == 0
              print("\r$(cl_i)/$(length(clients)) clients sifted")
          end
          client_tables[cl] = data[clients[cl], [:rssilog_rssi, :rssilog_antenna, :rssilog_timestamp2]]
      end
      println("")
      return client_tables
  end

  function read_client_data()
      c_path = "/Users/LTV/Documents/export_rssi_fixed.csv"

      println("Reading data")
      client = readtable(c_path)
      println("Process time")
      remove_milliseconds(client)
      println("Process IPs")
      ip_to_labels!(client)

      return client[:,[:rssilog_rssi, :rssilog_antenna, :rssilog_timestamp2, :session_sso]]
  end

  function save_clients(clients)
      path_to_write = "/Users/LTV/Work Drive/WiFi Positioning/rssi_data/town/data/clients_large_log/"

      for cl in keys(clients)
          path_to_write = "/Users/LTV/Work Drive/WiFi Positioning/rssi_data/town/data/clients_large_log/client_$(cl).csv"

          sort!(clients[cl], cols=[:rssilog_timestamp2])

          writetable(path_to_write, clients[cl])
      end
  end

end
