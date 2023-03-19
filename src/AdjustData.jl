export Reshaperaw!, OrderSlice!


# FUNCTION FOR REVERSING K-SPACE TRAJECTORIES TO ACCOUNT FOR BIPOLAR GRADIENTS
# might not work for Philips data
function Reshaperaw!(rawd::RawAcquisitionData)

    total_num = length(rawd.profiles)
    slices = zeros(Float32, total_num)

    for ii=1:total_num

        slices[ii] = rawd.profiles[ii].head.position[3]
        flags = digits(rawd.profiles[ii].head.flags, base = 2, pad=31)

        if flags[19] == true
            global noisemat = rawd.profiles[ii].data
            global noise_idx = ii
        end

        if flags[22] == true
            reverse!(rawd.profiles[ii].data, dims=1)
            rawd.profiles[ii].head.flags=rawd.profiles[ii].head.flags-(2^21)
        end
    
    end

    if @isdefined noise_idx
        deleteat!(rawd.profiles, noise_idx)
    end

    unique!(slices)
    slices_indx = sortperm(sortperm(slices))
    total_num = length(rawd.profiles)
    for ii = 1:total_num
        rawd.profiles[ii].head.idx.slice = slices_indx[rawd.profiles[ii].head.idx.slice+1]-1
    end

    if @isdefined noise_idx
        return slices_indx
    end

end

# FUNCTION TO SPATIALLY ORDER THE SLICES IN THE RAW DATA STRUCTURE
# might not work for Philips data

function OrderSlice!(rawd::RawAcquisitionData)

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