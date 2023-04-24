

function find_wrapped(nav, nav_time, trace, slices, TR)

    time = trace[:,1]
    trace_data = trace[:,2] ./ findmax(trace[:,2])[1] .* pi/2
    nav_norm = deepcopy(nav)
    nav_norm = nav_norm[1,1,:,:]

    nav_int = interpolate_mat(nav_norm, nav_time, trace_data, time, slices)


    # TEST TO BE REMOVED
    nav_slices = reshape(nav_norm, (slices * size(nav_norm,1)))
    nav_time_slices = reshape(nav_time, (slices * size(nav_norm,1)))
    nav_align = interpolate_mat(nav_slices, nav_time_slices, trace_data, time, 0)
    trace_time = align(nav_align, nav_time_slices, trace_data, time, TR)

    time_relevant = findall(x -> (x>(findmin(abs.(nav_time))[1]+TR) && x< (findmax(abs.(nav_time))[1]-TR)), time)
    correlation = cor(trace_data[time_relevant], nav_int[time_relevant,:])

    slices_index = findall(x -> x> 0.3, abs.(correlation))
    slices_index = [ x[2] for x in slices_index ]
    num_slices = size(slices_index,1)

    if num_slices > slices/3
        nav_slices = reshape(nav_norm[:,slices_index], (num_slices * size(nav_norm,1)))
        nav_time_slices = reshape(nav_time[1,1,:,slices_index], (num_slices * size(nav_norm,1)))
        nav_align = interpolate_mat(nav_slices, nav_time_slices, trace_data, time, 0)

        # align data and trace
        trace_time = align(nav_align, nav_time_slices, trace_data, time, TR)
    else
        trace_time = time
        nav_time_slices = nav_time
    end

    (nav_int, time, trace_data_red) = region_correlation(nav_int, time, trace_data, trace_time, nav_time_slices, TR)

    @mput nav_time trace_data trace_time slices
    mat"""
        for ii = 1:slices
            trace_data_inter(:,ii) = interp1(trace_time, trace_data, squeeze(nav_time(1,1,:,ii)));
        end
        %plot(squeeze(nav_time(1,1,:,1)), squeeze(trace_data_inter(:,1)), '.')
    """
    @mget trace_data_inter

    # renormalize trace and nav data
    for ii = 1:slices
        trace_data_inter[:,ii] = trace_data_inter[:,ii] ./ findmax(trace_data_inter[:,ii], dims =1)[1]
        nav_norm[:,ii] = nav_norm[:,ii] ./ findmax(nav_norm[:,ii], dims =1)[1]
    end

    # FLAG FIELD SHIFTS DIRECTION
    # compute navigator baseline
    nav_baseline = find_line(nav_norm, trace_data_inter, slices, "baseline")
    nav_topline = find_line(nav_norm, trace_data_inter, slices, "topline")
    correlation = ones(Int8, slices)

    for ii = 1:slices
        if nav_baseline[ii] <= nav_topline[ii]
            correlation[ii] = 1
        else
            correlation[ii] = -1
        end
    end

    # Invert slices with negative correlation to get a positive field shift
    for ii = 1:slices
        if correlation[ii] == -1
            nav_norm[:,ii] = -1 .* nav_norm[:,ii]
        end
    end

    # compute navigator baseline
    nav_baseline = find_baseline(nav_norm, trace_data_inter, slices)

    # reposition data to align the baseline
    for ii = 1:slices
        nav_norm[:,ii] = nav_norm[:,ii] .- nav_baseline[ii]
    end

    wrapped_points = find_wrapped_points(nav_norm, trace_data_red, trace_data_inter)

    mat"""
        clc
        clear
    """

    # return position wrapped points and field shift direction
    return wrapped_points, correlation

end



#FUNCTION TO INTERPOLATE THE DATA USING MATLAB
function interpolate_mat(nav_data, nav_time, trace_data, time, slices)
    
    # interpolate nav data to the trace time
    @mput nav_data nav_time trace_data time slices
    mat"""
    if not(slices == 0)
        for ii = 1:slices
            nav_int(:,ii) = interp1(squeeze(nav_time(1,1,:,ii)), squeeze(nav_data(:,ii)), time);
            %plot(time, nav_int, '.-');
        end
    
    else
        nav_int = interp1(squeeze(nav_time), squeeze(nav_data), time);
        %figure()
        %plot(time, nav_int);
        %hold on
        %plot(time, trace_data);
    end
    """
    @mget nav_int
    return nav_int

end



function align(nav_align, nav_time_slices, trace_data, time, TR)

    # align the two recordings
    TR = Int(TR)
    tmp = 1e10
    shift = 0
    time_relevant = findall(x -> (x>(findmin(abs.(nav_time_slices))[1]) && x< (findmax(abs.(nav_time_slices))[1])), time)
    for ii = -TR:1:TR
        resid_std = std(nav_align[time_relevant] - trace_data[time_relevant.+ii])
        if resid_std < tmp
            tmp=resid_std
            shift = ii
        end
    end

    print(string(shift, " ")) # test to check the time shift
    trace_time = time .- shift

    return trace_time

end


function find_baseline(nav_norm, trace_data_inter, slices)

    nav_baseline = zeros(Float64, slices)

    for ii = 1:slices
        line = (findmax(trace_data_inter[:,ii])[1] - findmin(trace_data_inter[:,ii])[1]) .*0.3 + findmin(trace_data_inter[:,ii])[1]
        relevant = findall(x -> x< line, trace_data_inter[:,ii])
        nav_norm_base = nav_norm[relevant,ii]
        nav_baseline[ii] = mean(nav_norm_base)
    end
    
    return nav_baseline

end


function find_line(nav_norm, trace_data_inter, slices, line_flag)

    nav_line = zeros(Float64, slices)

    for ii = 1:slices
        if line_flag == "baseline"
            line_base = -0.8
            line_top = 0
        elseif line_flag == "topline"
            line_base = 0
            line_top = 0.8
            
        end
        relevant = findall(x -> line_base <= x < line_top, trace_data_inter[:,ii])
        nav_norm_base = nav_norm[relevant,ii]
        nav_line[ii] = mean(nav_norm_base)
    end
    
    return nav_line

end


function region_correlation(nav_int, time, trace_data, trace_time, nav_time_slices, TR)

    # select relevant data region (only relevant for computing correlation)
    TR_safe = TR + TR*0.2
    time_relevant = findall(x -> (x>=(findmin(abs.(nav_time_slices))[1]+TR_safe) && x<=(findmax(abs.(nav_time_slices))[1]-TR_safe)), time)
    nav_int = nav_int[time_relevant,:] #interpolated nav relevant
    time = time[time_relevant] #interpolated nav time relevant
    time_relevant = findall(x -> (x>=(findmin(abs.(nav_time_slices))[1]+TR_safe) && x<=(findmax(abs.(nav_time_slices))[1]-TR_safe)), trace_time)
    #trace_time = trace_time[time_relevant]
    trace_data_red = trace_data[time_relevant] # trace data relevant
    if size(time,1) > size(trace_data_red, 1)
        time = time[1:size(trace_data_red, 1)]
        nav_int = nav_int[1:size(trace_data_red, 1),:]
    elseif size(time,1) < size(trace_data_red, 1)
        trace_data_red = trace_data_red[1:size(time,1)]
    end

    return nav_int, time, trace_data_red

end


function find_wrapped_points(nav_norm, trace_data_red, trace_data_inter)

    # check wrapped points
    wrapped_points = zeros(Int8, size(nav_norm))
    wrap_min = findmax(trace_data_red)[1] - ((findmax(trace_data_red)[1] - findmin(trace_data_red)[1]) .*0.36)
    idx_pos = findall(x -> x >= wrap_min, trace_data_inter)
    nav_add2pi = findall(x->x<0, nav_norm[idx_pos])
    wrapped_points[idx_pos[nav_add2pi]] .= 1

    return wrapped_points

end


function wrap_corr(nav, wrapped_points, correlation, slices)

    nav =  invert_nav(nav, correlation, slices)
    wrapped_points = reshape(wrapped_points, (1, 1, size(wrapped_points)...))
    idx_pos = findall(x->x==1, wrapped_points)
    nav[idx_pos] = nav[idx_pos] .+ (2*pi)
    nav =  invert_nav(nav, correlation, slices)

    return nav

end


function invert_nav(nav, correlation, slices)

    for ii = 1:slices
        if correlation[ii] == -1
            nav[:,:,:,ii] = -1 .* nav[:,:,:,ii]
        end
    end

    return nav

end