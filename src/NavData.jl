export additionalNavInput, navOutput

mutable struct additionalNavInput

    numslices::Int64
    numechoes::Int64
    numsamples::Int64
    numlines::Int64
    TE_nav::Float64
    dt_nav::Float64
    freq_enc_FoV::Union{Array{Int64}, Nothing}
    freq_enc_samples::Union{Array{Int64}, Nothing}
    phase_enc_samples::Union{Array{Int64}, Nothing}
    nav_time::Union{Array{Float64, 2}, Nothing}
    noisemat::Array{Complex{Float32}, 2}
    trace::Union{Matrix{Float64}, Nothing}
    centerline::Union{Vector{Float64}, Nothing}

end

"""
    Data = additionalNavInput(
        noisemat::Array{Complex{Float32}, 2},
        rawData::RawAcquisitionData,
        acqData::AcquisitionData,
        acqMap::Union{AcquisitionData, Nothing} = nothing,
        nav_time::Union{Array{Complex{Float32}, 2}, Nothing} = nothing,
        trace::Union{Matrix{Float64}, Nothing} = nothing,
        centerline::Union{Vector{Float64}, Nothing} = nothing)

Construct the additional data structure that is needed as imput to navCorr!

# Arguments
* `noisemat::Array{Complex{Float32}, 2}` - noise data obtained with ExtractNoiseData!
* `rawData::RawAcquisitionData` - raw data structure obtained loading raw data with MRIReco.jl
* `acqData::AcquisitionData` - acquisition data structure obtained converting raw data with MRIReco.jl

# Optional arguments with default value = nothing
* `acqMap::Union{AcquisitionData, Nothing} = nothing`       - acquisition data structure obtained converting reference data with MRIReco.jl
* `nav_time::Union{Array{Complex{Float32}, 2}, Nothing}`    - time stamps for the navigator data obtained with ExtractNavigator (in ms from the beginning of the day)
* `trace::Union{Matrix{Float64}, Nothing}`                  - respiratory trace time stamps and values in matrix with two colunms (1:time [ms], 2:trace).
                                                                Include time points before and after the image acquisition (at least 2 s).
* `centerline::Union{Vector{Float64}, Nothing}`             - coordinates of the spinal cord ceterline obtained with callSCT

MRIReco reference: https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.28792
"""
function additionalNavInput(
        noisemat::Array{Complex{Float32}, 2},
        rawData::RawAcquisitionData,
        acqData::AcquisitionData,
        acqMap::Union{AcquisitionData, Nothing} = nothing,
        nav_time::Union{Array{Float64, 2}, Nothing} = nothing,
        trace::Union{Matrix{Float64}, Nothing} = nothing,
        centerline::Union{Vector{Float64}, Nothing} = nothing    
    )

    numslices = numSlices(acqData)
    numechoes = numContrasts(acqData)
    numsamples = acqData.encodingSize[1]
    numlines = convert(Int64, size(acqData.kdata[1],1)/numsamples)

    ii=1
    while rawData.profiles[ii].head.user_int[8] < rawData.profiles[ii+1].head.user_int[8]
        ii=ii+1
    end

    # set up nav timing
    dt_nav = convert(Float64, rawData.profiles[ii-1].head.sample_time_us) .* 1e-6
    TE_nav = rawData.profiles[ii].head.user_int[8] .* 1e-6 # get TE nav

    if !isnothing(acqMap) && !isnothing(acqData)
        (freq_enc_FoV, freq_enc_samples, phase_enc_FoV, phase_enc_samples) = Find_scaling_sensit(acqMap, acqData)
    end

    return additionalNavInput(numslices, numechoes, numsamples, numlines, TE_nav, dt_nav,
                freq_enc_FoV, freq_enc_samples, phase_enc_samples, nav_time, noisemat, trace, centerline) 

end

mutable struct navOutput

    navigator::Array{Float64, 4}
    centerline::Union{Array{Float64, 1}, Nothing}
    correlation::Union{Array{Float64, 1}, Matrix{Float64}, Nothing}
    wrapped_points::Union{Array{Int8, 2}, Nothing}
    
end