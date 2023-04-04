export OrderSlices!, ExtractNoiseData!, ReverseBipolar!, RemoveRef!, CopyTE!, AdjustSubsampleIndices!

"""
    OrderSlices!(rawData::RawAcquisitionData)

Spatially order the slices in the raw data structure.

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

Extract the acquisition flags from raw data profiles.
Return a 31 elements vector for each profile.

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

Extract and return the noise acquisition from the raw data.
The noise acquisition is one of the profiles with slice = 0, contrast = 0, repetition = 0.

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

Reflect the raw data profiles for bipolar acquisition.

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

        if flags[ii,22] == true
            reverse!(rawData.profiles[ii].data, dims=1)
            rawData.profiles[ii].head.flags=rawData.profiles[ii].head.flags-(2^21)
        end

    end

end


"""
    RemoveRef!(rawData::RawAcquisitionData, slices::Union{Vector{Int64}, Nothing}, echoes::Union{Vector{Int64}, Nothing})

Remove reference data that are acquired with the phase stabilization on Siemens scanners.
Not solid to recalls.

# Arguments
* `rawData::RawAcquisitionData` - raw data structure obtained loading raw data with MRIReco.jl
* `slices::Union{Vector{Int64}, Nothing}` - vector containing the numbers of slices to be loaded with MRIReco.jl. Nothing loads all.
* `echoes::Union{Vector{Int64}, Nothing}` - vector containing the numbers of echoes to be loaded with MRIReco.jl. Nothing loads all.
"""
function RemoveRef!(rawData::RawAcquisitionData, slices::Union{Vector{Int64}, Nothing}, echoes::Union{Vector{Int64}, Nothing})

    numSlices = 0
    numEchoes = 0
    if slices === nothing
        numSlices = rawData.params["enc_lim_slice"].maximum+1
    else
        numSlices = size(slices, 1)
    end
    if echoes !== nothing
        if 0 in echoes
            numEchoes = size(echoes, 1) +1 # the navigator is saved as echo zero
        else
            numEchoes = size(echoes, 1)
        end
    else
        numEchoes = size(rawData.params["TE"],1)+1
    end

    #Apply this only if using phase stabilizaion
    removeIndx = numSlices*(numEchoes)
    deleteat!(rawData.profiles, 1:removeIndx)

end


"""
    CopyTE!(rawData::RawAcquisitionData, acqData::AcquisitionData)

Copy the TE values from the raw data structor to the acquisition data structor.

# Arguments
* `rawData::RawAcquisitionData` - raw data structure obtained loading raw data with MRIReco.jl
* `acqData::RawAcquisitionData` - acquisition data structure obtained converting raw data with MRIReco.jl
"""
function CopyTE!(rawData::RawAcquisitionData, acqData::AcquisitionData)

    for ii=1:size(acqData.kdata)[1]
        acqData.traj[ii].TE = rawData.params["TE"][ii]
    end

end


"""
    AdjustSubsampleIndices!(acqData::AcquisitionData)

Add subsamples indices in the acquisition data structure.
Needed when conveting data not acquired at the first repetition.

# Arguments
* `acqData::RawAcquisitionData` - acquisition data structure obtained converting raw data with MRIReco.jl
"""
function AdjustSubsampleIndices!(acqData::AcquisitionData)

    if isempty(acqData.subsampleIndices[1])
        for ii = 1:size(acqData.subsampleIndices)[1]
            acqData.subsampleIndices[ii]=1:size(acqData.kdata[1,1,1])[1]
        end
    end

end