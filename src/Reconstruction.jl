export Reconstruct, directreco

"""
    img = Reconstruct(acqd::AcquisitionData, sensit::Array{Complex{T},4}, noisemat::Union{Array{Complex{T}},Nothing} = nothing)

Call MRIReco.jl reconstruction function and return reconstructed image. Only single repetition in input.

# Arguments
* `acqData::RawAcquisitionData` - acquisition data structure obtained converting raw data with MRIReco.jl
* `sensit::Array{Complex{T},4}` - coil sensitivity map matric computed with CompSensit(acq::AcquisitionData, thresh = 0.135)
* `noisemat::Union{Array{Complex{T}},Nothing} = nothing` - noise data extracted from the raw datat structure with ExtractNoiseData!(rawData::RawAcquisitionData)

MRIReco reference: https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.28792
"""
function Reconstruct(acqd::AcquisitionData,
                    sensit::Array{Complex{T},4},
                    noisemat::Union{Array{Complex{T}},Nothing} = nothing) where {T} 

    params = Dict{Symbol, Any}()
    params[:reco] = "multiCoil"
    params[:solver] = "cgnr"
    params[:regularization] = "L2"
    params[:λ] = 1.e-2
    params[:iterations] = 10
    params[:reconSize] = (acqd.encodingSize[1],acqd.encodingSize[2])
    params[:estimateProfileCenter] = true
    params[:senseMaps] = sensit

    if !isnothing(noisemat)
        params[:noiseData] = noisemat
    end

    # Do reconstruction
    img = reconstruction(acqd, params)
    img = dropdims(img, dims = tuple(findall(size(img) .== 1)...))
    
    return img
end


"""
    img = directreco(acq::AcquisitionData)

Call MRIReco.jl reocnstruction function and return reconstructed image. Reconstruct coils separately.

# Arguments
* `acqData::RawAcquisitionData` - acquisition data structure obtained converting raw data with MRIReco.jl

MRIReco reference: https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.28792
"""
function directreco(acq::AcquisitionData)

    params = Dict{Symbol, Any}()
    params[:reco] = "direct"
    params[:reconSize] = (acq.encodingSize[1],acq.encodingSize[2])

    return reconstruction(acq, params)

end