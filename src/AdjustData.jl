export OrderSlices!, ExtractFlags, ExtractNoiseData!, ReverseBipolar!, RemoveRef!, CopyTE!, AdjustSubsampleIndices!

# FUNCTION TO SPATIALLY ORDER THE SLICES IN THE RAW DATA STRUCTURE
# might not work for Philips data

function OrderSlices!(rawd::RawAcquisitionData)

    total_num = length(rawd.profiles)
    slices = zeros(typeof(rawd.profiles[1].head.position[3]), total_num)

    for ii = 1:total_num
        slices[ii] = rawd.profiles[ii].head.position[3]
    end

    unique!(slices)
    slices_indx = sortperm(sortperm(slices))

    for ii = 1:total_num
        rawd.profiles[ii].head.idx.slice = slices_indx[rawd.profiles[ii].head.idx.slice+1]-1
    end

end


function ExtractFlags(rawd::RawAcquisitionData)

    total_num = length(rawd.profiles)
    flags = zeros(Int64, total_num, 31)

    for ii=1:total_num
        flags[ii,:] = digits(rawd.profiles[ii].head.flags, base = 2, pad=31)
    end

    return flags

end


function ExtractNoiseData!(rawData::RawAcquisitionData, flags::Array{Int64}) 

    total_num = length(rawData.profiles)
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

function ReverseBipolar!(rawd::RawAcquisitionData, flags::Array{Int64})

    total_num = length(rawd.profiles)

    for ii=1:total_num

        if flags[ii,22] == true
            reverse!(rawd.profiles[ii].data, dims=1)
            rawd.profiles[ii].head.flags=rawd.profiles[ii].head.flags-(2^21)
        end

    end

end


# FUNCTION FOR REMOVING THE REFERENCE DATA
function RemoveRef!(rawData::RawAcquisitionData, slices::Union{Int64, Nothing}, echoes::Union{Int64, Nothing})

    numSlices = 0
    numEchoes = 0
    if slices === nothing
        numSlices = rawData.params["enc_lim_slice"].maximum+1
    else
        slices == size(slices, 1)
    end
    if echoes !== nothing
        if 0 in echoes
            numEchoes = size(echoes, 1) +1
        else
            numEchoes = size(echoes, 1)
        end
    else
        numEchoes = size(rawData.params["TE"],1)+1
    end

    #Apply this only if using phase stabilizaion
    removeIndx = numSlices*(numEchoes) #*numberrep
    deleteat!(rawData.profiles, 1:removeIndx)

end


# FUNCTION TO COPY TE VALUES FROM RAW TO ACQ STRUCTURE
function CopyTE!(rawData::RawAcquisitionData, acqData::AcquisitionData)

    for ii=1:size(acqData.kdata)[1]
        acqData.traj[ii].TE = rawData.params["TE"][ii]
    end

end

function AdjustSubsampleIndices!(acqData::AcquisitionData)

    if isempty(acqData.subsampleIndices[1])
        for ii = 1:size(acqData.subsampleIndices)[1]
            acqd.subsampleIndices[ii]=1:size(acqd.kdata[1,1,1])[1]
        end
    end

end