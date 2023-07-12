function find_wrapped(nav::Array{Float64, 4}, nav_time::Array{Float64, 2}, trace::Array{Float64, 2}, slices::Int64, TR::Int64)

    time = trace[:,1]
    trace_data = trace[:,2] ./ findmax(trace[:,2])[1] .* pi/2
    trace_data = smooth_trace(time, trace_data)
    nav_norm = deepcopy(nav[1,1,:,:])
    nav_norm = smooth_nav(nav_time, nav_norm, slices)

    # Interpolate the trace values to the navigator time points for each slice
    trace_data_int = interpolate(trace_data, time, nav_time, slices)
    correlation = signalCorrelation(nav_norm, trace_data_int, slices, true)

    # Find navigator slices with higher correlationt to the trace
    max_corr = findmax(abs.(correlation))[1]
    std_corr = std(abs.(correlation))
    corr_relevant = findall(>(max_corr - 1.5*std_corr), abs.(correlation))
    size_corr_relevant = size(corr_relevant,1)

    # Invert navigator sign if the correlation is negative
    invertNavSign!(nav_norm, correlation, slices)

    # reshape in one vector the navigator signal from the slices with higher correlation, and align with the trace
    nav_align = reshape(nav_norm[:,corr_relevant], (size_corr_relevant * size(nav_norm,1)))
    nav_time_align = reshape(nav_time[:,corr_relevant], (size_corr_relevant * size(nav_time,1)))
    order = sortperm(nav_time_align)
    nav_align = nav_align[order]
    nav_time_align = nav_time_align[order]
    nav_align = smooth_trace(nav_time_align, nav_align)
    nav_align = interpolate(nav_align, nav_time_align, time)
    trace_time = align(nav_align, nav_time_align, trace_data, time, TR)

    # Invert navigator sign if the correlation is negative
    invertNavSign!(nav_norm, correlation, slices)

    # Interpolate the trace values to the navigator time points for each slice
    time_relevant = findall(x -> (x>(findmin(abs.(nav_time_align))[1]) && x< (findmax(abs.(nav_time_align))[1])), trace_time)
    trace_data_int = interpolate(trace_data[time_relevant], trace_time[time_relevant], nav_time, slices)

    # Compute correlation after alignemnt
    correlation = signalCorrelation(nav_norm, trace_data_int, slices, true)

    # how many filed changes?
    filter = digitalfilter(Lowpass(0.07, fs = 1), Butterworth(3))
    padd_size = size(correlation, 1)
    corr_padd = zeros(Float64, padd_size)
    corr_filt = cat(corr_padd .= correlation[1], correlation, corr_padd .= correlation[end], dims=1)
    corr_filt = filtfilt(filter, corr_filt)
    corr_filt = corr_filt[padd_size+1:end-padd_size,:]
    sign_corr_filt = sign.(corr_filt)
    field_change = 0
    for ii = 1:slices-1
        if sign_corr_filt[ii] != sign_corr_filt[ii+1]
            field_change = field_change +1
        end
    end

    corr_sign = sign.(correlation)

    if field_change == 0 || field_change == 2
        sign_corr = sign(mean(corr_sign))
        correlation = abs.(correlation) .* sign_corr
    elseif field_change == 1
        index_field_change = findlast(sign_corr_filt .== sign_corr_filt[1])[1]
        index_vector = collect(index_field_change - 2 : index_field_change + 2)
        filter!(x-> x != -1 && x != -2 && x != slices +1 && x != slices +2, index_vector)
        mean_val = 10
        index_smooth = 0
        index_corr = 0
        for ii in index_vector
            deviation = std(nav_norm[:,ii])
            meanval = mean(nav_norm[:,ii])
            remove_extreme = findall(x -> meanval - 0.8*deviation < x < meanval + 0.8*deviation, nav_norm[:,ii])
            mean_tmp = abs(mean(nav_norm[remove_extreme,ii]))
            if mean_tmp < mean_val
                mean_val = mean_tmp
                print(meanval)
                index_smooth = ii
            end
        end
        tmp = 0
        for ii in index_vector[1:end-1]
            if corr_sign[ii] != corr_sign[ii+1]
                if tmp == 0
                    index_corr = ii
                    tmp = 1
                end
            end
        end
        sign_corr_filt[1:index_corr] .= sign_corr_filt[1]
        sign_corr_filt[index_corr+1:end] .= sign_corr_filt[slices]
        correlation = abs.(correlation) .* sign_corr_filt
    end

    # Invert navigator sign if the correlation is negative
    invertNavSign!(nav_norm, correlation, slices)

    # find slices with possible wrapping
    possible_wrap_slices = ones(Bool, slices)
    for ii = 1:slices
        deviation = std(nav_norm[:,ii])
        meanval = mean(nav_norm[:,ii])
        remove_extreme = findall(x -> meanval - deviation < x < meanval + deviation, nav_norm[:,ii])
        max_interval = findmax(nav_norm[remove_extreme,ii])[1] - findmin(nav_norm[remove_extreme,ii])[1]
        if max_interval < 0.9
            possible_wrap_slices[ii] = false
        end
    end

    if mean(correlation) < 0.2
        possible_wrap_slices .== false
    end

    nowrap_slices = findall(possible_wrap_slices .== false)
    nav_norm[:, nowrap_slices] .= 1
    trace_data_int[:, nowrap_slices] .= 1

    # renormalize trace and nav data
    for ii = 1:slices
        trace_data_int[:,ii] = trace_data_int[:,ii] ./ findmax(trace_data_int[:,ii], dims =1)[1]
        nav_norm[:,ii] = nav_norm[:,ii] ./ findmax(nav_norm[:,ii], dims =1)[1]
    end

    # compute navigator baseline
    nav_baseline = find_baseline(nav_norm, trace_data_int, slices)

    # reposition data to align the baseline
    for ii = 1:slices
        nav_norm[:,ii] = nav_norm[:,ii] .- nav_baseline[ii]
    end

    time_relevant = findall(x -> (x>(findmin(abs.(nav_time_align))[1]) && x< (findmax(abs.(nav_time_align))[1])), trace_time)
    wrapped_points = find_wrapped_points(nav_norm, trace_data_int, trace_data[time_relevant])

    # return position wrapped points and field shift direction
    return wrapped_points, correlation

end

function smooth_trace(time::Array{Float64, 1}, trace_data::Array{Float64, 1})

    sampling_freq = size(time,1)/(time[end]-time[1])*1000
    filter = digitalfilter(Lowpass(0.7, fs = sampling_freq), Butterworth(3))
    return filtfilt(filter, trace_data)

end

function smooth_nav(nav_time::Array{Float64, 2}, nav_norm::Array{Float64, 2}, slices:: Int64)

    sampling_freq = size(nav_time, 1) * size(nav_time,2) / (findmax(abs.(nav_time))[1] - findmin(abs.(nav_time))[1]) *1000
    filter = digitalfilter(Highpass(0.5, fs = sampling_freq), Butterworth(3))
    padd_size = div(size(nav_norm, 1),4)
    nav_padd = zeros(Float64, padd_size, slices)
    nav_filt = cat(nav_padd, nav_norm, nav_padd, dims=1)
    for ii = 1:slices
        nav_filt[:,ii] = filtfilt(filter, nav_filt[:,ii])
    end
    nav_filt = nav_filt[padd_size+1:end-padd_size,:]
    return nav_filt

end


function interpolate(nav_norm::Union{Matrix{Float64}, Vector{Float64}},
                    nav_time::Union{Matrix{Float64}, Vector{Float64}},
                    time::Union{Matrix{Float64}, Vector{Float64}}, slices = 0)

    if ndims(nav_norm) == 1 && ndims(time) == 1
        nav_int = zeros(Float64, size(time,1))
        interp = DataInterpolations.LinearInterpolation(nav_norm, nav_time)
        nav_int = interp(time)
    elseif ndims(nav_norm) > 1 && ndims(time) == 1
        nav_int = zeros(Float64, size(time,1), size(nav_norm,2))
        for ii = 1:slices
            interp = DataInterpolations.LinearInterpolation(nav_norm[:,ii], nav_time[:,ii])
            nav_int[:,ii] = interp(time)
        end
    elseif ndims(nav_norm) == 1 && ndims(time) > 1
        nav_int = zeros(Float64, size(time))
        interp = DataInterpolations.LinearInterpolation(nav_norm, nav_time)
        for ii = 1:slices
            nav_int[:,ii] = interp(time[:,ii])
        end
    end

    return nav_int
end

function signalCorrelation(nav_int::Array{Float64, 2}, trace_data::Array{Float64, 2}, slices::Int64, allData::Bool)

    corr = ones(Float64, slices)
    for ii = 1:slices
        if allData == false
            deviation = std(trace_data[:,ii])
            meanval = mean(trace_data[:,ii])
            remove_extreme = findall(x -> x < meanval + deviation, trace_data[:,ii])
            corr[ii] = cor(trace_data[remove_extreme,ii], nav_int[remove_extreme,ii])
        else
            corr[ii] = cor(trace_data[:,ii], nav_int[:,ii])
        end
    end
    index_nan = isnan.(corr)
    index_nan = findall(index_nan .== 1)
    corr[index_nan] .= 0.1

    return corr
end

function invertNavSign!(nav::Union{Array{Float64, 2}, Array{Float64, 4}}, correlation::Union{Array{Float64, 1}, Matrix{Float64}}, slices::Int64)

    corr_sign = sign.(correlation)
    dimensions = ndims(nav)
    if dimensions == 2
        for ii = 1:slices
            nav[:,ii] = nav[:,ii] * corr_sign[ii]
        end
    elseif dimensions == 4
        for ii = 1:slices
            nav[:,:,:,ii] = nav[:,:,:,ii] * corr_sign[ii]
        end
    end
end

#FUNCTION TO ALIGN THE TRACE AND NAVIGATOR DATA
function align(nav_align::Array{Float64, 1}, nav_time_align::Array{Float64, 1}, trace_data::Array{Float64, 1}, time::Array{Float64, 1}, TR::Int64)

    time_relevant = findall(x -> (x>(findmin(abs.(nav_time_align))[1]) && x< (findmax(abs.(nav_time_align))[1])), time)
    delay = alignsignals(trace_data[time_relevant], nav_align[time_relevant])[2]
    trace_time = time
    if delay < TR/2
        print(delay)
        trace_time = circshift(time, delay)
    end

    return trace_time

end


function find_baseline(nav_norm::Array{Float64, 2}, trace_data_int::Array{Float64, 2}, slices::Int64)

    nav_baseline = zeros(Float64, slices)

    for ii = 1:slices
        deviation = std(trace_data_int[:,ii])
        if deviation == 0
            deviation = 1
        end
        meanval = mean(trace_data_int[:,ii])
        remove_extreme = findall(x -> meanval - deviation < x < meanval + deviation, trace_data_int[:,ii])
        trace_remove_extreme = trace_data_int[remove_extreme, ii]
        nav_remove_extreme = nav_norm[remove_extreme,ii]
        max = findmax(trace_remove_extreme)[1]
        min = findmin(trace_remove_extreme)[1]
        line = (max - min) .* 0.2 + min
        relevant = findall(x -> x< line, trace_remove_extreme)
        nav_baseline[ii] = mean(nav_remove_extreme[relevant])
    end
    
    return nav_baseline

end



#FUNCTION TO FIND THE WRAPPED POINTS
function find_wrapped_points(nav_norm::Array{Float64, 2}, trace_data_int::Array{Float64, 2}, trace_data_red::Array{Float64, 1})

    deviation = std(trace_data_red)
    meanval = mean(trace_data_red)
    remove_extreme = findall(x -> meanval - deviation < x < meanval + deviation, trace_data_red)
    wrapped_points = zeros(Int8, size(nav_norm))
    wrap_min = findmax(trace_data_red[remove_extreme])[1] - ((findmax(trace_data_red[remove_extreme])[1] - findmin(trace_data_red[remove_extreme])[1]) .*0.3)
    idx_pos = findall(x -> x >= wrap_min, trace_data_int)
    nav_add2pi = findall(x->x< -0.22, nav_norm[idx_pos])
    wrapped_points[idx_pos[nav_add2pi]] .= 1

    return wrapped_points

end