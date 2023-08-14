"""
    find_wrapped(nav::Array{Float64, 4}, nav_time::Array{Float64, 2}, trace::Array{Float64, 2}, slices::Int64, TR::Int64)

Identify the position of the wrapped points in the navigator phase estimates. The respiratory belt recording is necessary.
Return the position of the wrapped points and the correlation between each navigator slice and the trace data.

# Arguments
* `nav::Array{Float64, 4}` - navigator phase estimates
* `nav_time::Array{Float64, 2}` - navigator data time stamps in ms from the beginning of the day, for each slice
* `trace::Array{Float64, 2}` - physiological trace recording. Two columns vector. The first column contains the time stamps in ms from the beginning of the day
* `slices::Int64` - number of slices
* `TR::Int64` - acqusition repetition time (TR)
"""
function find_wrapped(nav::Array{Float64, 4}, nav_time::Array{Float64, 2}, trace::Array{Float64, 2}, slices::Int64, TR::Int64)

    time = trace[:,1]
    trace_data = trace[:,2] ./ findmax(trace[:,2])[1] .* pi/2
    trace_data = smooth_trace(time, trace_data)
    nav_norm = deepcopy(nav[1,1,:,:])
    nav_norm = smooth_nav(nav_time, nav_norm, slices)

    # Interpolate the trace values to the navigator time points for each slice
    trace_data_int = interpolate(trace_data, time, nav_time, slices)
    correlation = signalCorrelation(nav_norm, trace_data_int, slices, true)

    # Invert navigator sign if the correlation is negative
    invertNavSign!(nav_norm, correlation, slices)

    # Find navigator slices with higher correlation to the trace
    max_corr = findmax(abs.(correlation))[1]
    std_corr = std(abs.(correlation))
    if isnan(std_corr)
        std_corr = 1
    end
    corr_relevant = findall(>(max_corr - 1.5*std_corr), abs.(correlation))
    size_corr_relevant = size(corr_relevant,1)

    # reshape in one vector the navigator signal from the slices with higher correlation, and align with the trace
    nav_align = reshape(nav_norm[:,corr_relevant], (size_corr_relevant * size(nav_norm,1)))
    nav_time_align = reshape(nav_time[:,corr_relevant], (size_corr_relevant * size(nav_time,1)))
    order = sortperm(nav_time_align)
    nav_align = nav_align[order]
    nav_time_align = nav_time_align[order]
    nav_align = smooth_trace(nav_time_align, nav_align) # smooth the signal after combining multiple slices
    nav_align = interpolate(nav_align, nav_time_align, time)
    trace_time = align(nav_align, nav_time_align, trace_data, time, TR)

    # Invert navigator sign if the correlation is negative
    invertNavSign!(nav_norm, correlation, slices)

    # Interpolate the trace values to the navigator time points for each slice
    time_relevant = findall(x -> (x>(findmin(abs.(nav_time_align))[1]) && x< (findmax(abs.(nav_time_align))[1])), trace_time)
    trace_data_int = interpolate(trace_data[time_relevant], trace_time[time_relevant], nav_time, slices)

    # Compute correlation after alignemnt
    correlation = signalCorrelation(nav_norm, trace_data_int, slices, true)

    # allow for only one sign change in the correlation across slices
    # adjust the correlation sign consequently
    # this is effective only if the number of slices is higher than 5
    if slices > 5
        correlation = find_field_changes(correlation, slices)
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

    if mean(correlation) < 0.2 # consider the trace inaccurate
        possible_wrap_slices .== false
    end

    nowrap_slices = findall(possible_wrap_slices .== false)
    nav_norm[:, nowrap_slices] .= 1
    trace_data_int[:, nowrap_slices] .= 1

    # renormalize nav data
    for ii = 1:slices
        nav_norm[:,ii] = nav_norm[:,ii] ./ findmax(nav_norm[:,ii], dims =1)[1]
    end

    # compute navigator baseline
    nav_baseline = find_baseline(nav_norm, trace_data_int, slices)

    # reposition data to align the baseline
    for ii = 1:slices
        nav_norm[:,ii] = nav_norm[:,ii] .- nav_baseline[ii]
    end

    wrapped_points = find_wrapped_points(nav_norm, trace_data_int, slices)

    # return position wrapped points and field shift direction
    return wrapped_points, correlation

end

"""
    trace_data = smooth_trace(time::Array{Float64, 1}, trace_data::Array{Float64, 1})

Smooth the physiological trace recording using a butterworth low-pass filter (cut-off frequency 0.7Hz, 3 poles)

# Arguments
* `time::Array{Float64, 1}` - time in ms from the beginning of the day for the belt recording
* `trace_data::Array{Float64, 1}` - belt recording
"""
function smooth_trace(time::Array{Float64, 1}, trace_data::Array{Float64, 1})

    sampling_freq = size(time,1)/(time[end]-time[1])*1000
    filter = digitalfilter(Lowpass(0.7, fs = sampling_freq), Butterworth(3))
    
    return filtfilt(filter, trace_data)

end

"""
    nav_norm = smooth_nav(nav_time::Array{Float64, 2}, nav_norm::Array{Float64, 2}, slices:: Int64)

Remove the low frequencies components from the navigatior phase estimate using a butterworth high-pass filter (cut-off frequency 0.5Hz, 3 poles)

# Arguments
* `nav_time::Array{Float64, 2}` - navigator data time stamps in ms from the beginning of the day, for each slice
* `nav_norm::Array{Float64, 1}` - navigator phase estimates
* `slices::Int64` - number of slices
"""
function smooth_nav(nav_time::Array{Float64, 2}, nav_norm::Array{Float64, 2}, slices::Int64)

    sampling_freq = size(nav_time, 1) * size(nav_time,2) / (findmax(abs.(nav_time))[1] - findmin(abs.(nav_time))[1]) *1000
    filter = digitalfilter(Highpass(0.15, fs = sampling_freq), Butterworth(5))

    padd_size = div(size(nav_norm, 1),4)
    nav_padd = zeros(Float64, padd_size, slices)
    nav_filt = cat(nav_padd, nav_norm, nav_padd, dims=1)

    for ii = 1:slices
        nav_filt[:,ii] = filtfilt(filter, nav_filt[:,ii])
    end

    nav_filt = nav_filt[padd_size+1:end-padd_size,:]

    return nav_filt

end

"""
    inter_vector = interpolate(nav_norm::Union{Matrix{Float64}, Vector{Float64}},
                   nav_time::Union{Matrix{Float64}, Vector{Float64}},
                   time::Union{Matrix{Float64}, Vector{Float64}}, slices = 0)

Interpolate the first input vector with time stamps specified in the second input to the time points specified in the third input.
Return the interpolartion result.

# Arguments
* `nav_norm::Union{Matrix{Float64}, Vector{Float64}}` - vector or matrix to be interpolated
* `nav_time::Union{Matrix{Float64}, Vector{Float64}}` - time points of the original data
* `time::Union{Matrix{Float64}, Vector{Float64}}` - desired time points
* `slices::Int64` - number of slices
"""
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

"""
    corr = signalCorrelation(nav_int::Array{Float64, 2}, trace_data::Array{Float64, 2}, slices::Int64, allData = true)

Return the correlation between the respiratory belt recording and the navigator field variations estimates.
The vector must be interpolated to the same time points before calling this function.
Return correlation = 0.1 if data is Nan.

# Arguments
* `nav_int::Array{Float64, 2}` - navigator phase estimes
* `trace_data::Array{Float64, 2}` - physiological trace recording from the respiratory belt
* `slices::Int64` - number of slices
* `allData::Bool` - use all the time points if true. Use only the lower time points in the belt trace if false, hopefully escluding wrapped points in the navigator
"""
function signalCorrelation(nav_int::Array{Float64, 2}, trace_data::Array{Float64, 2}, slices::Int64, allData = true)

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


"""
    invertNavSign!(nav::Union{Array{Float64, 2}, Array{Float64, 4}}, correlation::Union{Array{Float64, 1}, Matrix{Float64}}, slices::Int64)

Invert the navigator phase estimates sign if the correlation between the respiratory trace and the navigator esimates is negative.

# Arguments
* `nav::Array{Float64, 2}` - navigator phase estimes
* `correlation::Union{Array{Float64, 1}, Matrix{Float64}}` - correlation vector between each navigator slice and the respiratory belt recording
* `slices::Int64` - number of slices
"""
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


"""
    trace_time = align(nav_align::Array{Float64, 1}, nav_time_align::Array{Float64, 1}, trace_data::Array{Float64, 1}, time::Array{Float64, 1}, TR::Int64)

Align the signal in the first imput (time stamps in the second imput) to the signal in the third imput (time stamps in the fourth input). acquisition TR in the last input.
Use the finddelay function from DSP.jl, find the peak of the signals cross-correlation.
Return the new time vector for the signal in the third input.

# Arguments
* `nav_align::Array{Float64, 2}` - navigator phase estimes reshaped in one vector
* `nav_time_align::Array{Float64, 1}` - time stamps for the navigator phase estimates in ms from the beginning of the day
* `trace_data::Array{Float64, 1}` - respiratory belt recording  in ms from the beginning of the day
* `time::Array{Float64, 1}` - time stamps for the respiratory belt recording in se
* `TR::Int64` - acquisition repetition time
"""
function align(nav_align::Array{Float64, 1}, nav_time_align::Array{Float64, 1}, trace_data::Array{Float64, 1}, time::Array{Float64, 1}, TR::Int64)

    time_relevant = findall(x -> (x>(findmin(abs.(nav_time_align))[1]) && x< (findmax(abs.(nav_time_align))[1])), time)
    delay = alignsignals(trace_data[time_relevant], nav_align[time_relevant])[2]
    trace_time = time

    if  -TR/2 < delay < TR/2
        trace_time = circshift(time, delay)
    end

    return trace_time

end


"""
    correlation = find_field_changes(correlation::Union{Array{Float64, 1}, Matrix{Float64}})

Inhale air can lead to both positive and negative field variations depensing by the vertebral level.
There are two regions where the field variations change sign, at the lungs extremities.
It is reasonable to assume that MRI using a commercial spinal coil can not allow to record both these regions in the same acquisition.
Therefore, only one field change in the correlation sign across slices should be allowed.
This function works only if the number of slices is bigger than 5.

# Arguments
* `correlation::Union{Array{Float64, 1}, Matrix{Float64}}` - correlation vector between each navigator slice and the respiratory belt recording
* `slices::Int64` - number of slices
"""
function find_field_changes(correlation::Union{Array{Float64, 1}, Matrix{Float64}}, slices::Int64)

    # allow for only one change in the sign on the correlation values across slices
    filter = digitalfilter(Lowpass(0.07, fs = 1), Butterworth(3))
    padd_size = size(correlation, 1)
    corr_padd = zeros(Float64, padd_size)
    corr_filt = cat(corr_padd .= correlation[1], correlation, corr_padd .= correlation[end], dims=1)

    # Wrap the filtering op in try-catch because filtfilt requires a minimum signal length
    try
        corr_filt = filtfilt(filter, corr_filt)
    catch e
        if isa(e,BoundsError)
            println("Please choose more slices, not enough to check field change across slices")
        end
    end

    # Count the field changes (sign changes)
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
        
        # init counters
        index_corr = 0
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

    return correlation

end


"""
    nav_baseline = find_baseline(nav_norm::Array{Float64, 2}, trace_data_int::Array{Float64, 2}, slices::Int64)

Find navigator baseline corresponding to the end of the exhalation. Return the baseline coordinate.

# Arguments
* `nav_norm::Array{Float64, 2}` - navigator phase estimes
* `trace_data_int::Array{Float64, 2}` - trace data smoothed, aligned and interpolated to the navigator time points for each slice.
* `slices::Int64` - number of slices
"""
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


"""
    wrapped_points = find_wrapped_points(nav_norm::Array{Float64, 2}, trace_data_int::Array{Float64, 2}, slices::Int64)

Find wrapped points comparing the breathing related oscillations measured with the respiratory belt and the navigator readout.
Return a binary array, with the same size as nav_norm and 1 if the point is idenfied as wrapped.

# Arguments
* `nav_norm::Array{Float64, 2}` - navigator phase estimes
* `trace_data_int::Array{Float64, 2}` - trace data smoothed, and interpolated to the navigator time points for each slice
* `slices::Int64` - number of slices
"""
function find_wrapped_points(nav_norm::Array{Float64, 2}, trace_data_int::Array{Float64, 2}, slices::Int64)

    wrapped_points = zeros(Int8, size(nav_norm))

    for ii = 1:slices

        deviation = std(trace_data_int[:,ii])

        if deviation == 0
            deviation = 1
        end

        meanval = mean(trace_data_int[:,ii])
        remove_extreme = findall(x -> meanval - deviation < x < meanval + deviation, trace_data_int[:,ii])
        wrap_min = findmax(trace_data_int[remove_extreme,ii])[1] - ((findmax(trace_data_int[remove_extreme,ii])[1] - findmin(trace_data_int[remove_extreme,ii])[1]) .*0.28)
        idx_pos = findall(x -> x >= wrap_min, trace_data_int[:,ii])
        nav_add2pi = findall(x->x< -0.2, nav_norm[idx_pos,ii])
        wrapped_points[idx_pos[nav_add2pi],ii] .= 1

    end

    return wrapped_points

end