export OrderSlices!, ExtractNoiseData!, ReverseBipolar!, RemoveRef!, CopyTE!, AdjustSubsampleIndices!, ExtractNavigator, ExtractFlags, selectEcho!, selectSlice!

"""
    OrderSlices!(rawData::RawAcquisitionData)

Spatially order the slices in the MRIReco.jl raw data structure.
The slices are ordered basing on the position coordinates saved in each profile.
If these are not present the slices can not be ordered.

MRIReco reference: https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.28792

# Arguments
* `rawData::RawAcquisitionData` - raw data structure obtained loading raw data with MRIReco.jl
"""
function OrderSlices!(rawData::RawAcquisitionData)

    total_num = length(rawData.profiles)
    slices = zeros(typeof(rawData.profiles[1].head.position[3]), total_num)

    for ii = 1:total_num
        slices[ii] = rawData.profiles[ii].head.position[3]
    end

    unique!(slices)
    slices_indx = sortperm(sortperm(slices))

    for ii = 1:total_num
        index = rawData.profiles[ii].head.position[3] .== slices
        rawData.profiles[ii].head.idx.slice = slices_indx[index][1]-1
    end

end


"""
    flags = ExtractFlags(rawData::RawAcquisitionData) 

Extract the acquisition flags from the MRIReco.jl raw data profiles.
Return a 31-element vector for each profile.

MRIReco reference: https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.28792

# Arguments
* `rawData::RawAcquisitionData` - raw data structure obtained loading raw data with MRIReco.jl
"""
function ExtractFlags(rawData::RawAcquisitionData)

    total_num = length(rawData.profiles)
    flags = zeros(Int64, total_num, 31)

    for ii=1:total_num
        flags[ii,:] = digits(rawData.profiles[ii].head.flags, base = 2, pad=31)
    end

    return flags

end


"""
    noisemat = ExtractNoiseData!(rawData::RawAcquisitionData, flags::Array{Int64})

Extract and return the noise acquisition from the MRIReco.jl raw data.
The noise acquisition is usually the first profile with slice = 0, contrast = 0, repetition = 0.
The noise profile should have the 19th flag element equal to 1. Check this with ExtractFlags if errors occur.

MRIReco reference: https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.28792

# Arguments
* `rawData::RawAcquisitionData` - raw data structure obtained loading raw data with MRIReco.jl
"""
function ExtractNoiseData!(rawData::RawAcquisitionData)

    flags = ExtractFlags(rawData)
    total_num = length(rawData.profiles)
    
    if total_num != size(flags, 1)
        @error "size of flags and number of profiles in rawData do not match"
    end

    noisemat = Matrix{typeof(rawData.profiles[1].data)}

    for ii=1:total_num

        # noise acquisition flag in ISMRMRD format
        if flags[ii,19] == true
            noisemat = rawData.profiles[ii].data
            deleteat!(rawData.profiles, ii)
            break
        end

    end

    return noisemat

end


"""
    ReverseBipolar!(rawData::RawAcquisitionData)

Reflect the MRIReco.jl raw data profiles for bipolar acquisition.
    
MRIReco reference: https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.28792

# Arguments
* `rawData::RawAcquisitionData` - raw data structure obtained loading raw data with MRIReco.jl
"""
function ReverseBipolar!(rawData::RawAcquisitionData)

    flags = ExtractFlags(rawData)
    total_num = length(rawData.profiles)

    if total_num != size(flags, 1)
        @error "size of flags and number of profiles in rawData do not match"
    end

    for ii=1:total_num

        # reflect line flag in ISMRMRD format
        if flags[ii,22] == true
            reverse!(rawData.profiles[ii].data, dims=1)
            rawData.profiles[ii].head.flags=rawData.profiles[ii].head.flags-(2^21)
        end

    end

end


"""
    RemoveRef!(rawData::RawAcquisitionData, slices::Union{Vector{Int64}, Nothing}, echoes::Union{Vector{Int64}, Nothing})

Remove reference data that are not useful for the navigator-based correction from acquisitions with phase stabilization on Siemens scanners.
Make sure that this is needed on your data checking the time stamps with mapVBVD in Matlab.
Not robust to repeated calls, modifies rawData.

MRIReco reference: https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.28792
mapVBVD reference: https://github.com/CIC-methods/FID-A/blob/master/inputOutput/mapVBVD/README.md

# Arguments
* `rawData::RawAcquisitionData` - raw data structure obtained loading raw data with MRIReco.jl
"""
function RemoveRef!(rawData::RawAcquisitionData)

    numSlices = rawData.params["enc_lim_slice"].maximum + 1
    numEchoes = size(rawData.params["TE"],1) + 1

    #Apply this only if using phase stabilizaion
    removeIndx = numSlices*(numEchoes)
    deleteat!(rawData.profiles, 1:removeIndx)

end


"""
    CopyTE!(rawData::RawAcquisitionData, acqData::AcquisitionData)

Copy the TE values from the MRIReco.jl raw data structure to the acquisition data structure.

MRIReco reference: https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.28792

# Arguments
* `rawData::RawAcquisitionData` - raw data structure obtained loading raw data with MRIReco.jl
* `acqData::AcquisitionData` - acquisition data structure obtained converting raw data with MRIReco.jl
"""
function CopyTE!(rawData::RawAcquisitionData, acqData::AcquisitionData)

    for ii=1:size(acqData.kdata)[1]
        acqData.traj[ii].TE = deepcopy(rawData.params["TE"][ii])
    end

end


"""
    AdjustSubsampleIndices!(acqData::AcquisitionData)

Add subsamples indices in the MRIReco.jl acquisition data structure.
Needed when converting data not acquired in the first repetition.

MRIReco reference: https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.28792

# Arguments
* `acqData::AcquisitionData` - acquisition data structure obtained converting raw data with MRIReco.jl
"""
function AdjustSubsampleIndices!(acqData::AcquisitionData)

    if isempty(acqData.subsampleIndices[1])

        for ii = 1:size(acqData.subsampleIndices)[1]
            acqData.subsampleIndices[ii]=1:size(acqData.kdata[1,1,1])[1]
        end

    end

end


"""
    (nav, nav_time) = ExtractNavigator(rawData::RawAcquisitionData, slices::Union{Vector{Int64}, Nothing})

Extract the navigator profiles from the MRIReco.jl raw data structure.
These are registered with the same indices (contract, slice, encoding step) as the image data for the first echo time.
Return a navigator array and a navigator time array. The navigator array has four dimensions in the following order: k-space samples, coils, k-space lines, slices.
Effective only if the navigator profile was acquired after the first image profile.
    
MRIReco reference: https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.28792

# Arguments
* `rawData::RawAcquisitionData` - raw data structure obtained loading raw data with MRIReco.jl
"""
function ExtractNavigator(rawData::RawAcquisitionData)

    @info "The navigaotr extraction is effective only if the navigator profile was acquired after the first image profile."

    total_num = length(rawData.profiles)
    numberslices = rawData.params["enc_lim_slice"].maximum +1
    contrasts = zeros(Int64, total_num)
    slices = zeros(Int64, total_num)
    lines = zeros(Int64, total_num)

    for ii = 1:total_num

        contrasts[ii] = rawData.profiles[ii].head.idx.contrast
        slices[ii] = rawData.profiles[ii].head.idx.slice
        lines[ii] = rawData.profiles[ii].head.idx.kspace_encode_step_1

    end

    # keep only the indexes of data saved in the first echo (this includes navigator)
    contrastsIndx = findall(x->x==0, contrasts)
    slices = slices[contrastsIndx]
    lines = lines[contrastsIndx]

    nav = zeros(typeof(rawData.profiles[1].data[1,1]), size(rawData.profiles[1].data)[1], size(rawData.profiles[1].data)[2], rawData.params["reconSize"][2], numberslices)

    nav_time = zeros(Float64, rawData.params["reconSize"][2], numberslices)

    #Odd indexes are data first echo, Even indexes are navigator data
    for ii = 2:2:length(slices)

        nav[:,:,lines[ii]+1,slices[ii]+1] = rawData.profiles[contrastsIndx[ii]].data
        nav_time[lines[ii]+1,slices[ii]+1] = rawData.profiles[contrastsIndx[ii]].head.acquisition_time_stamp

    end

    #Remove the rows filled with zeroes
    lines = unique(lines) .+1
    nav = nav[:,:,lines,:]
    nav_time = nav_time[lines,:]

    return nav, nav_time
    #navigator[k-space samples, coils, k-space lines, slices]

end

"""
    SelectEcho!(acqd, idx_echo)

Extract one or more echoes from the acquisition data structure

# Arguments
* `acqd::AcquisitionData` - acquisition data structure obtained converting raw data with MRIReco.jl
* `idx_echo::Vector{Int64}` - vector containing the indexes of the echoes to be selected (starting from 0)

MRIReco reference: https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.28792
"""
function selectEcho!(acqd::AcquisitionData, idx_echo::Vector{Int64})

    if !isempty(idx_echo)

        contrasts = size(acqd.kdata)[1]
        indices = Vector{Int64}(undef, contrasts)

        for ii=1:contrasts
            indices[ii] = ii
        end

        deleteat!(indices, idx_echo)
        deleteat!(acqd.subsampleIndices, indices)
        deleteat!(acqd.traj, indices)
        acqd.kdata = acqd.kdata[idx_echo,:,:]

    end
end

"""
    SelectSlice!(acqd, nav, nav_time, idx_slice)

Extract one or more echoes from the acquisition data structure

# Arguments
* `acqd::AcquisitionData` - acquisition data structure obtained converting raw data with MRIReco.jl
* `idx_slice::Vector{Int64}` - vector containing the indexes of the slices to be selected (starting from 0, downer slice)

# Optional arguments with default value = nothing
* `nav::Union{Array{Complex{T}, 4}, Nothin} = nothing` - navigator profiles obtained with the ExtractNavigator function
* `nav_time::Union{Array{Complex{Float32}, 2}, Nothing}` - time stamps for the navigator data obtained with ExtractNavigator (in ms from the beginning of the day)

MRIReco reference: https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.28792
"""
function selectSlice!(acqd::AcquisitionData, idx_slice::Vector{Int64}, nav::Union{Array{Complex{T}, 4}, Nothing} = nothing, nav_time::Union{Array{Float64, 2}, Nothing} = nothing) where {T}

    # get kdata from slice
    acqd.kdata = acqd.kdata[:,idx_slice,:]

    if !isnothing(nav) && !isnothing(nav_time)

        nav = nav[:,:,:,idx_slice]
        nav_time = nav_time[:,idx_slice]

    end

end