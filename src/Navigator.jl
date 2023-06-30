export NavCorr!, comp_weights, comp_centerline, TE_corr, apply_corr!, remove_ref_ph!

function NavCorr!(nav::Array{Complex{T}, 4}, acqData::AcquisitionData, params::Dict{Symbol, Any}, addData::additionalNavInput) where{T}
    
    #navigator[k-space samples, coils, k-space lines, slices]
    # compute the navigator fourier transform in the readout direction
    centerline = nothing
    if params[:corr_type] != "knav"
        nav = ifftshift(ifft(fftshift(nav, [1]), [1]), [1])
        #noisemat = fftshift(fft(ifftshift(noisemat, [1]), [1]), [1])
        nav_center = div(addData.numsamples, 2) # change name, center is a function
        if params[:use_SCT] == true
            centerline = comp_centerline(addData)
            for ii = 1:addData.numslices
                nav[:,:,:,ii] = circshift(nav[:,:,:,ii], nav_center-centerline[ii])
            end
        end
        buff = round(Int64, params[:FFT_interval] / acqData.fov[1] * addData.numsamples / 2)
        nav = nav[nav_center-buff+1:nav_center+buff,:,:,:]
    end

    remove_ref_ph!(nav, addData.numlines, 1) # remove the reference phase
    navabs = abs.(nav)
    noisestd = std(addData.noisemat, dims=[1]).^2

    # Compute weights for the coils average
    weights = comp_weights(navabs, noisestd, addData.numlines, addData.numslices)
    nav = sum(weights .* nav, dims=(1,2,)) # coils and lines average

    cartes_index = findall(x -> isnan(x), nav) # remove this, ugly fix subject 01 rep 2
    phMean = nav./abs.(nav)
    phMean[cartes_index] .= 0 # remove this, ugly fix subject 01 rep 2
    phMean = angle.(mean(phMean, dims=(3,))) # Compute mean phase
    nav = nav./exp.(im*phMean) # Recenter phase of time series
    nav = angle.(nav) # compute navigator phase
    nav[cartes_index] .= 0 # remove this, ugly fix subject 01 rep 2

    correlation = nothing
    wrapped_points = nothing
    
    if params[:corr_type] == "FFT_unwrap"
        (wrapped_points, correlation) = find_wrapped(nav, addData.nav_time, addData.trace, addData.numslices, addData.TR)
        nav = wrap_corr(nav, wrapped_points, correlation, addData.numslices)
    end

    nav_return = deepcopy(nav)
    
    # Correct for different TEs
    nav = TE_corr(nav, acqData, addData.dt_nav, addData.TE_nav, addData.numsamples, addData.numechoes)
    nav = exp.(im*nav)
    # Apply the correction to the data
    apply_corr!(nav, acqData, addData.numechoes,addData.numlines, addData.numsamples, addData.numslices)
    return navOutput(nav_return, centerline, correlation, wrapped_points)

end

function comp_weights(navabs, noisestd, lines::Int64, slices::Int64)

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


function comp_centerline(addData::additionalNavInput)

    freq_enc_ref_res = addData.freq_enc_FoV[1] / addData.freq_enc_samples[1]
    freq_enc_img_res = addData.freq_enc_FoV[2] / addData.freq_enc_samples[2]
    freq_enc_FoV_disc = Int64((addData.freq_enc_FoV[1] - addData.freq_enc_FoV[2]) / freq_enc_ref_res / 2)
    start_voxel = div(addData.freq_enc_samples[1] - addData.phase_enc_samples[1], 2)
    centerline = round.(addData.centerline)
    centerline = centerline .+ start_voxel
    centerline = centerline .- freq_enc_FoV_disc
    centerline = centerline .* freq_enc_ref_res ./ freq_enc_img_res
    centerline = floor.(centerline)

    return centerline

end

function TE_corr(nav::Array{T, 4}, acqd::AcquisitionData, dt_nav::Float64, TE_nav::Float64, numsamples::Int64, numechoes::Int64) where {T}

    nav = nav ./ TE_nav
    nav = repeat(nav, outer=(numsamples,numechoes,1,1))
    t_nav = ([1:numsamples;] .- (numsamples/2+1)) .* dt_nav
    for ii=1:numechoes
        for ll=1:numsamples
            nav[ll,ii,:,:] = nav[ll,ii,:,:] .* (acqd.traj[ii].TE .* 1e-3 + t_nav[ll])
        end
    end
    return nav

end

function apply_corr!(nav::Array{T, 4}, acqd::AcquisitionData, numechoes::Int64, numlines::Int64, numsamples::Int64, numslices::Int64) where {T}

    for ii = 1:numechoes
        for jj = 1:numlines
            for ll = 1:numsamples
                for mm = 1:numslices
                    acqd.kdata[ii,mm,1][(jj-1)*numsamples+ll,:] =
                     acqd.kdata[ii,mm,1][(jj-1)*numsamples+ll,:] ./ nav[ll,ii,jj,mm]
                end
            end
        end
    end

end

function remove_ref_ph!(nav::Array{Complex{T}, 4}, lines::Int64, index::Int64) where {T}

    phRef = exp.(im*angle.(nav[:,:,index,:])) # compute reference phase
    for ii=1:lines
        nav[:,:,ii,:] = nav[:,:,ii,:] ./ phRef # subctract reference phase
    end

end


function wrap_corr(nav::Array{Float64, 4}, wrapped_points::Array{Int8, 2}, correlation::Array{Float64, 1}, slices::Int64)

    invertNavSign!(nav, correlation, slices)
    wrapped_points = reshape(wrapped_points, (1, 1, size(wrapped_points)...))
    idx_pos = findall(x->x==1, wrapped_points)
    nav[idx_pos] = nav[idx_pos] .+ (2*pi)
    invertNavSign!(nav, correlation, slices)

    return nav
end