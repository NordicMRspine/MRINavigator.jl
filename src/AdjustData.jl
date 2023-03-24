export OrderSlices!

# FUNCTION TO SPATIALLY ORDER THE SLICES IN THE RAW DATA STRUCTURE
# might not work for Philips data

function OrderSlices!(rawd::RawAcquisitionData)

    total_num = length(rawd.profiles)
    slices = zeros(Float32, total_num)

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


function ExtractNoiseData!(rawd::RawAcquisitionData, flags::Array{Int64}) 

    total_num = length(rawd.profiles)
    noisemat = Matrix{ComplexF32}

    for ii=1:total_num

        if flags[ii,19] == true
            noisemat = rawd.profiles[ii].data
            deleteat!(rawd.profiles, ii)
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
function RemoveRef!(rawd, slices, echoes)

    #Apply this only if using phase stabilizaion
    removeIndx = slices*(echoes+1) #*numberrep
    deleteat!(rawd.profiles, 1:removeIndx)

end