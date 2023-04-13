export NavCorr!

# FUNCTION TO APPLY THE NAVIGATOR CORRECTION WITH FFT
function NavCorr!(nav::Array{Complex{T}, 4}, acqData::AcquisitionData, params::Dict{Symbol, Any}, addData::additionalDataStruct) where{T}
    
    #navigator[k-space samples, coils, k-space lines, slices]
    # compute the navigator fourier transform in the readout direction
    if params[:corr_type] != "knav"
        nav = ifftshift(ifft(fftshift(nav, [1]), [1]), [1])
        #noisemat = fftshift(fft(ifftshift(noisemat, [1]), [1]), [1])
        center = div(addData.numsamples, 2) # change name, center is a function
        if params[:use_SCT] == true
            centerline = comp_centerline(addData)
            for ii = 1:addData.numslices
                nav[:,:,:,ii] = circshift(nav[:,:,:,ii], center-centerline[ii])
            end
        end
        buff = Int64(params[:FFT_interval] / acqData.fov[1] * addData.numsamples / 2)
        nav = nav[center-buff+1:center+buff,:,:,:]
    end

    remove_ref_ph!(nav, addData.numlines, 1) # remove the reference phase
    navabs = abs.(nav)
    noisestd = std(addData.noisemat, dims=[1]).^2

    # Compute weights for the coils average
    weights = comp_weights(navabs, noisestd, addData.numlines, addData.numslices)
    nav = sum(weights .* nav, dims=(1,2,)) # coils and lines average

    phMean = nav./abs.(nav)
    phMean = angle.(mean(phMean, dims=(3,))) # Compute mean phase
    nav = nav./exp.(im*phMean) # Recenter phase of time series
    nav = angle.(nav) # compute navigator phase

    if corr_type == "FFT_narrow_wrap"
        (wrapped_points, correlation) = find_wrapped(nav, nav_time, trace, slices, TR)
        nav = wrap_corr(nav, wrapped_points, correlation, slices)
    end

    nav_return = deepcopy(nav)
    
    # Correct for different TEs
    nav = TE_corr(nav, acqd, dt_nav, TE_nav, samples, contrasts)
    nav = exp.(im*nav)
    # Apply the correction to the data
    apply_corr!(nav, acqd, contrasts, lines, samples, slices)
    return nav_return, centerline

end

function comp_weights(navabs, noisestd, lines, slices)

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


function comp_centerline(addData)

    freq_enc_ref_res = addData.freq_enc_FoV[1] / addData.freq_enc_samples[1]
    freq_enc_img_res = addData.freq_enc_FoV[2] / addData.freq_enc_samples[2]
    freq_enc_FoV_disc = Int64((addData.freq_enc_FoV[1] - addData.freq_enc_FoV[2]) / freq_enc_ref_res / 2)
    centerline = round.(addData.centerline)
    centerline = centerline .- freq_enc_FoV_disc
    centerline = centerline .* freq_enc_ref_res ./ freq_enc_img_res
    centerline = floor.(centerline)

    return centerline

end

# FUNCTION TO CORRECT FOR DIFFERENT TEs
function TE_corr(nav, acqd, dt_nav, TE_nav, samples, contrasts)

    nav = nav ./ TE_nav
    nav = repeat(nav, outer=(samples,contrasts,1,1))
    t_nav = ([1:samples;] .- (samples/2+1)) .* dt_nav
    for  ii=1:contrasts
        for ll=1:samples
            nav[ll,ii,:,:] = nav[ll,ii,:,:] .* (acqd.traj[ii].TE .* 1e-3 + t_nav[ll])
        end
    end
    return nav

end

function apply_corr!(nav, acqd, contrasts, lines, samples, slices)

    for ii = 1:contrasts
        for jj = 1:lines
            for ll = 1:samples
                for mm = 1:slices
                    acqd.kdata[ii,mm,1][(jj-1)*samples+ll,:] =
                     acqd.kdata[ii,mm,1][(jj-1)*samples+ll,:] ./ nav[ll,ii,jj,mm]
                end
            end
        end
    end

end

function remove_ref_ph!(nav, lines, index)

    phRef = exp.(im*angle.(nav[:,:,index,:])) # compute reference phase
    for ii=1:lines
        nav[:,:,ii,:] = nav[:,:,ii,:] ./ phRef # subctract reference phase
    end

end