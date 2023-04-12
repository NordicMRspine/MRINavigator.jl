mutable struct additionalDataStruct
    numslices::Int64
    numechoes::Int64
    numsamples::Int64
    numlines::Int64
    TR::Union{Int64, Nothing}
    TE_nav::Int64
    dt_nav::Float64
    freq_enc_FoV::Union{Int64, Nothing}
    freq_enc_samples::Union{Int64, Nothing}
    nav_time::Union{Array{Complex{Float32}, 2}, Nothing}
    noisemat::Array{Complex{Float32}, 2}
    trace::Union{Matrix{Float64}, Nothing}
    centerline::Union{Vector{Float64}, Nothing}
end


function additionalDataStruct(
    acqMap::AcquisitionData,
    acqData::AcquisitionData,
    nav_time::Array{Complex{Float32}, 2},
    noisemat::Array{Complex{Float32}, 2},
    trace::Matrix{Float64},
    centerline::Vector{Float64},
    rawData::RawAcquisitionData,
    ) where {T}

    numslices = size(acqData.kdata)[2]
    numechoes = size(acqData.kdata)[1]
    numsamples = acqData.encodingSize[1]
    numlines = convert(Int64, size(acqData.kdata[1],1)/samples)
    TR = rawData.params["TR"]
    ii=1
    while rawData.profiles[ii].head.user_int[8] < rawData.profiles[ii+1].head.user_int[8]
        ii=ii+1
    end
    dt_nav = convert(Float64, rawData.profiles[ii-1].head.sample_time_us) .* 1e-6
    TE_nav = rawd.profiles[ii].head.user_int[8] .* 1e-6 # get TE nav
    (freq_enc_FoV, freq_enc_samples) = Find_scaling_sensit(acqMap, acqData)


return Data = (numslices, numechoes, numsamples, numlines, TR, TE_nav, dt_nav,
                freq_enc_FoV, freq_enc_samples, nav_time, noisemat, trace, centerline)

end