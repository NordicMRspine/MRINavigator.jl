

function find_wrapped(nav::Array{Float32,4}, nav_time::Array{Float32, 2}, trace::Array{Float64, 1}, slices::Int64, TR::Int64) where {T}

    time = trace[:,1]
    trace_data = trace[:,2] ./ findmax(trace[:,2])[1] .* pi/2
    nav_norm = deepcopy(nav[1,1,:,:])

    # Interpolate the navigator values from each slice to trace sample points
    nav_int = interpolate(nav_norm, nav_time, time, slices)

    # find navigato time limits for each slice
    time_lim = zeros(Int64, slices, 2)
    for ii = 1:slices
        time_lim[ii,1] = findfirst(>=(nav_time[1,ii]), time)
        time_lim[ii,2] = findfirst(>=(nav_time[end,ii]), time)
    end

    # compute correlation between navgator signal and trace signal
    correlation = signalCorrelation(nav_int, trace_data, time_lim)
    max_corr = findmax(abs.(correlation))[1]
    std_corr = std(abs.(correlation))
    # Find slices with signal correlation within 1.5 sigma from the maximum
    corr_relevant = findall(>(max_corr - 1.5 * std_corr), abs.(correlation))
    size_corr_relevant = size(corr_relevant,1)
    # Invert navigator sign if the correlation is negative
    invertNavSign!(nav_norm, correlation, slices)

    # reshape the navigator signal from the slices with higher correlation in a vector and align with the trace
    nav_align = reshape(nav_norm[:,corr_relevant], (size_corr_relevant * size(nav_norm,1)))
    nav_time_align = reshape(nav_time[:,corr_relevant], (size_corr_relevant * size(nav_time,1)))
    order = sortperm(nav_time_align)
    nav_align = nav_align[order]
    nav_time_align = nav_time_align[order]
    nav_align = interpolate(nav_align, nav_time_align, time, 0) # zero flag for no slices
    trace_time = align(nav_align, nav_time, trace_data, time, TR)

    # Re-Invert navigator sign if the correlation is negative
    invertNavSign!(nav_norm, correlation, slices)
    # Compute correlation after alignemnt
    correlation = signalCorrelation(nav_int, trace_data, time_lim)




end



function interpolate(nav_norm::Array{Complex{T}, 2}, nav_time::Array{Float32, 2}, time::Array{Float64, 1}, slices::Int64)

    if slices == 0
        nav_int = zeros(Float64, size(time,1))
        interp = DataInterpolations.LinearInterpolation(nav_norm, nav_time)
        nav_int = interp(time)
    else
        nav_int = zeros(Float64, size(time,1), size(nav_norm,2))
        for ii = 1:slices
            interp = DataInterpolations.LinearInterpolation(nav_norm[:,ii], nav_time[:,ii])
            nav_int[:,ii] = interp(time)
        end
    end
    return nav_int
end

function signalCorrelation(nav_int::Array{Complex{T}, 2}, trace_data::Array{Float64, 1}, time_lim::Array{Int64, 2}) where {T}

    corr = ones(Float64, slices)
    for ii = 1:slices
        time_relevant = time_lim[ii,1] : time_lim[ii,2]
        corr[ii] = cor(trace_data[time_relevant], nav_int[time_relevant, ii])
    end

    return corr
end

function invertNavSign!(nav::Array{Complex{T}, 2}, correlation::Array{Float64, 1}, slices::Int64) where {T}

    corr_sign = sign.(correlation)

    for ii = 1:slices
        nav[:,ii] = nav[:,ii] * corr_sign[ii]
    end

end

#FUNCTION TO ALIGN THE TRACE AND NAVIGATOR DATA
function align(nav_align::Array{Complex{T}, 1}, nav_time_slices::Array{Float32, 1}, trace_data::Array{Float64, 1}, time::Array{Float64, 1}, TR::Int64)

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

    #print(string(shift, " ")) # test to check the time shift
    trace_time = time .- shift

    return trace_time

end