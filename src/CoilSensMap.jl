export CompSensit, ResizeSensit!, CompRoughMask


# FUNCTION TO COMPUTE THE COIL SENSITIVITY MAP
function CompSensit(acq::AcquisitionData, thresh = 0.138)

    sensit = espirit(acq,(6,6),30,eigThresh_1=0.005, eigThresh_2=0)
    slices = numSlices(acq)
    coils = numChannels(acq)
    # compute mask
    mask = CompRoughMask(acq, slices, thresh)

    for ii = 1:slices

        mask_slice = mask[:,:,ii]
        findConnectedComponent!(mask_slice)
        removeBehindBack!(mask_slice)
        homogeneousMask!(mask_slice)
        mask[:,:,ii] = mask_slice

    end
    for ii = 1:coils
        sensit[:,:,:,ii] = sensit[:,:,:,ii] .* mask
    end

    return sensit

end

function CompRoughMask(acq::AcquisitionData, slices::Int64, thresh = 0.138)

    img = directreco(acq)
    I_sum = sqrt.(sum(abs.(img) .^ 2, dims = 5)) .+ eps()
    I_sum = dropdims(I_sum, dims = tuple(findall(size(I_sum) .== 1)...))
    I_max = ones(Float64, slices)
    mask = zeros(size(I_sum))
    for ii = 1:slices
        I_max[ii] = maximum(abs.(I_sum[:,:,ii]))
        mask[findall(x -> x > thresh * I_max[ii], I_sum[:,:,ii]), ii] .= 1
    end

    return mask
end


function findConnectedComponent!(mask_slice::Array{T,2}) where {T}

    components = label_components(mask_slice)
    measured_area = component_lengths(components)
    measured_area = measured_area[2:end] #remove background component
    blob = findmax(measured_area)[2]
    cartes_index_blob = findall(x -> x!=blob, components)
    mask_slice[cartes_index_blob] .= 0

end

function removeBehindBack!(mask_slice::Array{T,2}) where{T}

    # remove noisy voxels on the left of the image
    # corresponding to the back of the subject
    density = sum(mask_slice, dims = 1)[1,:]
    # compute dervative of points density
    lines = div(size(mask_slice, 2),2)
    dder = zeros(Int64, lines)
    for ii = 1:lines
        dder[ii] = density[ii+1] - density[ii]
    end
    # put to zero everything behind the patient back
    position = findmax(dder)[2] -3
    for jj=1:position
        mask_slice[:,jj].=0
    end

end

function homogeneousMask!(mask_slice::Array{T,2}) where{T}

    cartes_index_slice = CartesianIndices(mask_slice)
    Bimask_slice = convert(BitMatrix, mask_slice)
    hull = convexhull(Bimask_slice)
    push!(hull, hull[1])
    inside = [inpolygon(p, hull; in=true, on=true, out=false) for p in cartes_index_slice]
    mask_slice[inside] .= 1

end

# FUNCTION TO INTERPOLATE AND REMOVE OUTFLINE OF THE COIL SENSITIVITY MAP
function ResizeSensit!(sensit::Array{Complex{T},4}, acqMap::AcquisitionData, acqData::AcquisitionData) where {T}

    # Define the relevant sensit region assuming the same slices center between ref and image data
    (freq_enc_FoV, freq_enc_samples, phase_enc_FoV, phase_enc_samples) = Find_scaling_sensit(acqMap, acqData)

    freq_enc_FoV_disc = Int64((freq_enc_FoV[1] - freq_enc_FoV[2]) / (freq_enc_FoV[1]/freq_enc_samples[1]) / 2)
    phase_enc_FoV_disc = Int64((phase_enc_FoV[1] - phase_enc_FoV[2]) / (phase_enc_FoV[1]/phase_enc_samples[1]) / 2)

    # Remove the sensit data that are not included in the image data FoV
    sensit = sensit[freq_enc_FoV_disc+1:end-freq_enc_FoV_disc, phase_enc_FoV_disc+1:end-phase_enc_FoV_disc, :, :]
    
    # Interpolate the sensit data to the image data
    cartes_index = findall(x -> x!=0, sensit)
    mask = zeros(Float32, size(sensit)) # compute the mask
    for ii in cartes_index
        mask[ii] = 1
    end
    # Linear interpolation
    sensit = mapslices(x ->imresize(x, (freq_enc_samples[2], phase_enc_samples[2])), sensit, dims=[1,2])
    mask = mapslices(x ->imresize(x, (freq_enc_samples[2], phase_enc_samples[2])), mask, dims=[1,2])
    # Remove interpolated outline
    cartes_index = findall(x -> x!=1, mask)
    for ii in cartes_index
        mask[ii] = 0
    end
    sensit = mask .* sensit

end

# FUNCTION TO INTERPOLATE AND REMOVE OUTFLINE OF THE COIL SENSITIVITY MAP
function Find_scaling_sensit(acqMap::AcquisitionData, acqData::AcquisitionData) where {T}

    freq_enc_FoV = [acqMap.fov[1], acqData.fov[1]]
    freq_enc_samples = [acqMap.encodingSize[1], acqData.encodingSize[1]]
    phase_enc_FoV = [acqMap.fov[2], acqData.fov[2]]
    phase_enc_samples = [acqMap.encodingSize[2], acqData.encodingSize[2]]

    return freq_enc_FoV, freq_enc_samples, phase_enc_FoV, phase_enc_samples

end