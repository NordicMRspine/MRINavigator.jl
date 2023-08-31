export NavCorr!, comp_centerline_pos, wrap_corr!, TE_corr!, apply_corr!

"""
    navOutput = NavCorr!(nav::Array{Complex{T}, 4}, acqData::AcquisitionData, params::Dict{Symbol, Any}, addData::additionalNavInput) where {T}

Compute the navigator-based correction and apply it to the acquisition data. Multiple pipelines are available: "knav", "FFT" and "FFT_unwrap".
Return navigator trace, spinal cord centerline in the reconstructed image coordinates, 
Correlation between navigator and belt data for each slice and position of wrapped points for each slices.
Please choose the pipeline using the corr_type filed in the params dictionary.

# Arguments
* `nav::Array{Complex{T}, 4}` - navigator profiles obtained with the ExtractNavigator function
* `acqData::AcquisitionData` - acquisition data structure obtained converting raw data with MRIReco.jl
* `params::Dict{Symbol, Any}` - navigator correction paramerters dictionary
* `addData::additionalNavInput` - mandatory additional data structure obtained with the constructor: additionalNavInput

MRIReco reference: https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.28792

"""
function NavCorr!(nav::Array{Complex{T}, 4}, acqData::AcquisitionData, params::Dict{Symbol, Any}, addData::additionalNavInput) where{T}
    
    # navigator[k-space samples, coils, k-space lines, slices]
    # compute the navigator fourier transform in the readout direction, only for FFT case
    centerline = nothing
    corr_type = split(params[:corr_type], "_")
    if corr_type[1] != "knav"
        nav = ifftshift(ifft(fftshift(nav, [1]), [1]), [1])
        # noisemat = fftshift(fft(ifftshift(noisemat, [1]), [1]), [1])

        nav_center = div(addData.numsamples, 2)
        if params[:use_centerline] == true
            centerline = comp_centerline_pos(addData)
            for ii = 1:addData.numslices
                nav[:,:,:,ii] = circshift(nav[:,:,:,ii], nav_center-centerline[ii])
            end
        end

        buff = round(Int64, params[:FFT_interval] / acqData.fov[1] * addData.numsamples / 2)

        # neglect points outside the interval of interest
        nav = nav[nav_center-buff+1:nav_center+buff,:,:,:]

    end

    remove_ref_ph!(nav, addData.numlines, 1) # remove the reference phase
    noisestd = std(addData.noisemat, dims=[1]).^2

    # Compute weights for the coils average
    weights = comp_weights(abs.(nav), noisestd, addData.numlines, addData.numslices)
    nav = sum(weights .* nav, dims=(1,2,)) # coils and samples average for each line

    # Compute navigator phase
    cartes_index = findall(x -> isnan(x), nav)
    nav[cartes_index] .= eps()
    phMean = nav./abs.(nav)
    phMean = angle.(mean(phMean, dims=(3,))) # Compute mean phase over the lines
    nav = nav./exp.(im*phMean) # Recenter phase of time series
    nav = angle.(nav) # compute navigator phase

    correlation = nothing
    wrapped_points = nothing
    
    corr_type = split(params[:corr_type], "_")
    if size(corr_type, 1) == 2
        if corr_type[2] == "unwrap"
            (wrapped_points, correlation) = find_wrapped(nav, addData.nav_time, addData.trace, addData.numslices)
            nav = wrap_corr!(nav, wrapped_points, correlation, addData.numslices)
        end
    end

    nav_return = deepcopy(nav)
    
    # Correct for different TEs
    nav = TE_corr!(nav, acqData, addData.dt_nav, addData.TE_nav, addData.numsamples, addData.numechoes)
    nav = exp.(im*nav)

    # Apply the correction to the data
    apply_corr!(nav, acqData, addData.numechoes,addData.numlines, addData.numsamples, addData.numslices)

    return navOutput(nav_return, centerline, correlation, wrapped_points)

end


"""
    weights = comp_weights(navabs::Array{T, 4}, noisestd::Matrix{T}, lines::Int64, slices::Int64) where {T}

Compute and return the weights for coils and sample average for each slice. These should be used to extract one value for each navigator profile.

# Arguments
* `navabs::Array{T, 4}` - module of the navigator data
* `noisestd::Matrix{T}` - standard deviation of the noise acquisition for each coil
* `lines::Int64` - number of lines or data profiles
* `slices::Int64` - number of slices
"""
function comp_weights(navabs::Array{T, 4}, noisestd::Matrix{T}, lines::Int64, slices::Int64) where {T}

    # weights[points, coils, lines, slices]
    coils = size(navabs, 2)
    weights = zeros(size(navabs))
    for ii=1:coils
        weights[:,ii,:,:] = navabs[:,ii,:,:] ./ noisestd[1,ii]
    end

    weightsnorm = sum(weights, dims=(1,2))

    for ii=1:lines
        for ll=1:slices
            weights[:,:,ii,ll] = weights[:,:,ii,ll] ./ weightsnorm[1,1,ii,ll]
        end
    end
    return weights

end

"""
    centerline = comp_centerline_pos(addData::additionalNavInput)

Convert and return centerline position from the reference data cordinate to the acquisition data coordinates (number of voxels).

# Arguments
* `addData::additionalNavInput` - mandatory additional data structure obtained with the constructor: additionalNavInput
"""
function comp_centerline_pos(addData::additionalNavInput)

    # Compute resolution and disc
    freq_enc_ref_res = addData.freq_enc_FoV[1] / addData.freq_enc_samples[1]
    freq_enc_img_res = addData.freq_enc_FoV[2] / addData.freq_enc_samples[2]
    freq_enc_FoV_disc = Int64((addData.freq_enc_FoV[1] - addData.freq_enc_FoV[2]) / freq_enc_ref_res / 2)

    start_voxel = div(addData.freq_enc_samples[1] - addData.phase_enc_samples[1], 2)
    
    # Compute centerline
    centerline = round.(addData.centerline)
    centerline = centerline .+ start_voxel
    centerline = centerline .- freq_enc_FoV_disc
    centerline = centerline .* freq_enc_ref_res ./ freq_enc_img_res
    centerline = floor.(centerline)

    return centerline

end

"""
    nav = TE_corr!(nav::Array{T, 4}, acqd::AcquisitionData, dt_nav::Float64, TE_nav::Float64, numsamples::Int64, numechoes::Int64) where {T}

Compute the phase value for the navigator correction basing on the exact acquisition time of each data sample in the line and for each echo.
Return a four-dimensional navigator array.

# Arguments
* `nav::Array{T, 4}` - phase estimates obtained from the navigator data
* `acqData::AcquisitionData` - acquisition data structure obtained converting raw data with MRIReco.jl
* `dt_nav::Float64` - time interval between two samples in the frequency encoding direction
* `TE_nav::Float64` - echo time of the navigator readout
* `numsamples::Int64` - number of samples for each profile
* `numechoes::Int64` - number of echoes

MRIReco reference: https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.28792
"""
function TE_corr!(nav::Array{T, 4}, acqd::AcquisitionData, dt_nav::Float64, TE_nav::Float64, numsamples::Int64, numechoes::Int64) where {T}

    # Set up navigator phase timing
    nav = nav ./ TE_nav
    nav = repeat(nav, outer=(numsamples,numechoes,1,1))
    t_nav = ([1:numsamples;] .- (numsamples/2+1)) .* dt_nav
    
    # Compute navigator phase
    for ii=1:numechoes
        for ll=1:numsamples
            nav[ll,ii,:,:] = nav[ll,ii,:,:] .* (acqd.traj[ii].TE .* 1e-3 + t_nav[ll])
        end
    end
    
    return nav

end

"""
    apply_corr!(nav::Array{T, 4}, acqd::AcquisitionData, numechoes::Int64, numlines::Int64, numsamples::Int64, numslices::Int64) where {T}

Apply the navigator-based correction to the acquisition data structure obtained loading the raw data with MRIReco.jl.
After applying the correction the image should be reconstructed. Use the reconstruct function.

# Arguments
* `nav::Array{T, 4}` - phase estimates obtained from the navigator data
* `acqd::AcquisitionData` - acquisition data structure obtained converting raw data with MRIReco.jl
* `numechoes::Int64` - number of echoes
* `numlines::Int64` - number of lines (profiles) for each slice and echo
* `numsamples::Int64` - number of samples for each profile
* `numslices::Int64` - number of slices

MRIReco reference: https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.28792
"""
function apply_corr!(nav::Array{T, 4}, acqd::AcquisitionData, numechoes::Int64, numlines::Int64, numsamples::Int64, numslices::Int64) where {T}

    # Loop over dims of the acquisition
    for ii = 1:numechoes
        for jj = 1:numlines
            for ll = 1:numsamples
                for mm = 1:numslices

                    # Apply correction
                    acqd.kdata[ii,mm,1][(jj-1)*numsamples+ll,:] = acqd.kdata[ii,mm,1][(jj-1)*numsamples+ll,:] ./ nav[ll,ii,jj,mm]

                end
            end
        end
    end

end

"""
    remove_ref_ph!(nav::Array{Complex{T}, 4}, lines::Int64, index::Int64) where {T}

Subtract the phase of a navigator profile to all the other profiles in the slice. The input index defines which is the reference profile.

# Arguments
* `nav::Array{Complex{T}, 4}` - navigator profiles
* `lines::Int64` - lines or profiles number
* `index::Int64` - index of the reference profile to be subtracted

"""
function remove_ref_ph!(nav::Array{Complex{T}, 4}, lines::Int64, index::Int64) where {T}
    
    # compute reference phase
    phRef = exp.(im*angle.(nav[:,:,index,:]))

    # subtract reference phase for each line
    for ii = 1:lines
        nav[:,:,ii,:] = nav[:,:,ii,:] ./ phRef 
    end

end

"""
    wrap_corr!(nav::Array{Float64, 4}, wrapped_points::Array{Int8, 2}, correlation::Union{Array{Float64, 1}, Matrix{Float64}}, slices::Int64)

Unwrap the wrapped points identified with the find_wrapped funtion. These functions can be used only if physiological recording is available.

# Arguments
* `nav::Array{T, 4}` - phase estimates obtained from the navigator data
* `wrapped_points::Array{Int8, 2}` - position of the wrapped points, output of find_wrapped
* `correlation::Union{Array{Float64, 1}` - correlation values between the physiological recording the navigator estimates for each slice. Output of find_wrapped
* `slices::Int64` - number of slices

"""
function wrap_corr!(nav::Array{Float64, 4}, wrapped_points::Array{Int8, 2}, correlation::Union{Array{Float64, 1}, Matrix{Float64}}, slices::Int64)

    invertNavSign!(nav, correlation, slices)
    wrapped_points_local = reshape(wrapped_points, (1, 1, size(wrapped_points)...))
    idx_pos = findall(x->x==1, wrapped_points_local)
    nav[idx_pos] = nav[idx_pos] .+ (2*pi)
    invertNavSign!(nav, correlation, slices)

    return nav
end