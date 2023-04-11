mutable struct NavDataStruct

    rawData::RawAcquisitionData,
    acqMap::AcquisitionData,
    acqData::AcquisitionData,
    nav::Array{Complex{Float32}, 4},
    nav_time::Array{Complex{Float32}, 2},
    noisemat::Array{Complex{Float32}, 2},
    trace::Matrix{Float64},
    centerline::Vector{Float64}

    function NavDataStruct(rawData::RawAcquisitionData,
                      acqMap::AcquisitionData,
                      acqData::AcquisitionData,
                      nav::Array{Complex{Float32}, 4},
                      nav_time::Array{Complex{Float32}, 2},
                      noisemat::Array{Complex{Float32}, 2},
                      trace::Matrix{Float64},
                      centerline::Vector{Float64},
                      ) where {T}

        return new(rawData, acqMap, acqData, nav, nav_time, noisemat, trace, centerline)

    end
end