export NavCorr!

# FUNCTION TO APPLY THE NAVIGATOR CORRECTION WITH FFT
function NavCorr!(nav::Array{Complex{T}, 4}, acqData::AcquisitionData, params::Dict{Symbol, Any}, addData::additionalDataStruct) where{T}

    # Get the initial parameters
    #navigator[k-space samples, coils, k-space lines, slices]
    (coils, slices, contrasts, samples, lines, TR, dt_nav, TE_nav) = get_init_params(acqd, rawd)

    # compute the navigator fourier transform in the readout direction
    if corr_type != "knav" && corr_type != "none"

        nav = ifftshift(ifft(fftshift(nav, [1]), [1]), [1])

        #noisemat = fftshift(fft(ifftshift(noisemat, [1]), [1]), [1])
        center = div(samples, 2) # change name, center is a function
        if config.use_SCT == true
            centerline = comp_centerline(acq, acqd, centerline)
            for ii = 1:slices
                nav[:,:,:,ii] = circshift(nav[:,:,:,ii], center-centerline[ii])
            end
        end
        buff = Int64(config.FFT_interval / acqd.fov[1] * samples / 2)
        nav = nav[center-buff+1:center+buff,:,:,:]
        
    end

    remove_ref_ph!(nav,lines,1) # remove the reference phase

    navabs = abs.(nav)
    noisestd = std(noisemat, dims=[1]).^2

    # Compute weights for the coils average
    weights = comp_weights(navabs, noisestd, lines, slices)
    
    nav = sum(weights .* nav, dims=(1,2,)) # coils and lines average
    #nav = sum(nav, dims=(1,2,))

    cartes_index = findall(x -> isnan(x), nav) # remove this, ugly fix subject 01 rep 2
    phMean = nav./abs.(nav)
    phMean[cartes_index] .= 0 # remove this, ugly fix subject 01 rep 2
    phMean = angle.(mean(phMean, dims=(3,))) # Compute mean phase
    nav = nav./exp.(im*phMean) # Recenter phase of time series
    nav = angle.(nav) # compute navigator phase

    # Check the wrapping
    #nav = wrap_corr_laura(nav, slices)
    #nav = wrap_corr(nav, wrap_points1, 1)
    #nav = wrap_corr(nav, wrap_points0, 0)

    nav[cartes_index] .= 0 # remove this, ugly fix subject 01 rep 2

    if corr_type == "FFT_narrow_wrap"
        (wrapped_points, correlation) = find_wrapped(nav, nav_time, trace, slices, TR)
        nav = wrap_corr(nav, wrapped_points, correlation, slices)

        path = config.path_res
        lable = config.lable
        @mput wrapped_points
        @mput path
        @mput lable
    
        mat"""
            figure()
            t = tiledlayout(6,3);
            for ii=1:18
                nexttile
                plot(1:204, squeeze(wrapped_points(:,ii)), '.-')
                xlim([-2,206])
                title(strcat("slice", string(ii)))
            end
            saveas(gcf, append(path,'/wrapped_', lable, '.fig'))
        """ 
    end

    nav_return = deepcopy(nav)
    
    # Correct for different TEs
    nav = TE_corr(nav, acqd, dt_nav, TE_nav, samples, contrasts)

    #nav = repeat(nav, outer=(1,1,1,samples))
    nav = exp.(im*nav)
    
    # Apply the correction to the data
    apply_corr!(nav, acqd, contrasts, lines, samples, slices)
    return nav_return, centerline

end

# FUNCTION TO GET THE INITIAL DATA PARAMETERS TO RUN THE NAVIGATOR CORRECTION
function get_init_params(acqd, rawd)

    coils = size(acqd.kdata[1,1,1])[2]
    slices = size(acqd.kdata)[2]
    contrasts = size(acqd.kdata)[1]
    samples = acqd.encodingSize[1]
    lines = convert(Int64, size(acqd.kdata[1],1)/samples)
    TR = rawd.params["TR"]
    ii=1
    while rawd.profiles[ii].head.user_int[8] < rawd.profiles[ii+1].head.user_int[8]
        ii=ii+1
    end
    dt_nav = convert(Float64, rawd.profiles[ii-1].head.sample_time_us) .* 1e-6
    TE_nav = rawd.profiles[ii].head.user_int[8] .* 1e-6 # get TE nav

    return coils, slices, contrasts, samples, lines, TR, dt_nav, TE_nav
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


function comp_centerline(acq, acqd, centerline)

    # Consider only navigator and noise data in the spinal cord region
    (freq_enc_FoV, freq_enc_samples) = Find_scaling_sensit(acq, acqd)
    freq_enc_ref_res = freq_enc_FoV[1] / freq_enc_samples[1]
    freq_enc_img_res = freq_enc_FoV[2] / freq_enc_samples[2]
    freq_enc_FoV_disc = Int64((freq_enc_FoV[1] - freq_enc_FoV[2]) / freq_enc_ref_res / 2)
    #recon_size = raw.params["reconSize"][1]
    #freq_enc_samples_offset = div(freq_enc_samples[1] - recon_size, 2)
    centerline = round.(centerline)
    #centerline = centerline .+ freq_enc_samples_offset .-1
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


# FUNCTION TO APPLY THE CORRECTION TO K-SPACE DATA
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