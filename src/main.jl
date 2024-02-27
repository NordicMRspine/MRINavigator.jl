export runNavPipeline, saveNoise, loadRawData, convertRawToAcq

"""
    runNavPipeline(params::Dict{Symbol, Any})

Run the navigator pipeline. Return reconstructed image and navigator correction output (check NavCorr!).

# Arguments
* `params::Dict{Symbol, Any}` - MRINavigator parameter structure, check defaultNavParams() for info
"""
function runNavPipeline(params::Dict{Symbol, Any})
    
    findCenterline(params)
    saveNoise(params[:path_imgData], params[:path_noise]::String)
    rawData = loadRawData(params)

    @info "load noise"
    # load noise nacquisition
    noisemat = FileIO.load(params[:path_noise], "noisemat")

    @info "Extract navigator data. The time stamps are accurate only for Siemens data."
    @info "The navigator extraction is effective only if the navigator profile was acquired after the first image profile."
    (nav, nav_time) = ExtractNavigator(rawData)
    nav_time = nav_time .* 2.5 # seconds from beginning of the day (Siemens data only)

    acqData = convertRawToAcq(rawData)

    # slice and echo selection on acquisition data
    selectEcho!(acqData, params[:echoes])
    (nav, nav_time) = selectSlice!(acqData, params[:slices], nav, nav_time)

    @info "read ref data"
    # read reference data
    rawMap = RawAcquisitionData(ISMRMRDFile(params[:path_refData]))
    OrderSlices!(rawMap)
    acqMap = AcquisitionData(rawMap, estimateProfileCenter=true)

    @info "sensemaps"
    ## compute or load the coil sensitivity map
    if params[:comp_sensit]
        CompResizeSaveSensit(acqMap, acqData, params[:path_sensit], params[:mask_thresh])
    end

    #Load coil sensitivity
    sensit = FileIO.load(params[:path_sensit], "sensit")
    if !isnothing(params[:slices])
        sensit = reshape(sensit[:,:,params[:slices],:],(size(sensit,1), size(sensit,2),
            size(params[:slices],1), size(sensit,4)))
    end

    # Load centerline (ON LINUX: file is centerline.csv, ON WINDOWS AND MAC: is centerline.nii.csv)
    centerline = nothing
    if params[:use_centerline] == true
        try
            run(`cat /etc/os-release`, wait = true)
        catch e
            if isa(e, ProcessFailedException)
                centerline = CSV.read(params[:path_centerline] * "centerline.nii.csv", DataFrame, header=false)
            else
                centerline = CSV.read(params[:path_centerline] * "centerline.csv", DataFrame, header=false)
            end
        end
        centerline = centerline.Column1
        if !isnothing(params[:slices])
            centerline = centerline[params[:slices]]
        end
    end

    #Load trace
    trace = nothing
    if params[:corr_type] == "FFT_unwrap"
        trace = read(matopen(params[:path_physio]), "data")
    end

    @info "nav corr"
    # Navigator correction
    if params[:corr_type] != "none"
        addData = additionalNavInput(noisemat, rawData, acqData, acqMap, nav_time, trace, centerline)
        output = NavCorr!(nav, acqData, params, addData)
    end

    @info "recon"
    ## Reconstruct the data
    img = Reconstruct(acqData, sensit, noisemat)

    return output, img

end


"""
    saveNoise(path_imgData::String, path_noise::String)

Extract the noise acquisition form the image data and save it.  
Call ExtractNoiseData!, check this function for more info.

# Arguments
* `path_imgData::String` - path to the ISMRMRD file containing the image data
* `path_noise::String` - path where the noise file will be saved

ISMRMRD reference: https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.26089
"""
function saveNoise(path_imgData::String, path_noise::String)

    @info "Load first rep, save noise acquisition"
    # load the first repetition, slice and echo and save the noise acquisition for optimal results
    # the noise acquisition is saved in the first repetition only
    rawData = RawAcquisitionData(ISMRMRDFile(path_imgData),
            slice = 0, contrast = 0, repetition = 0)
    noisemat = ExtractNoiseData!(rawData)
    FileIO.save(path_noise,"noisemat",noisemat)

end


"""
    loadRawData(params::Dict{Symbol, Any})

Load the raw data file saved in ISMRMRD format in julia using MRIReco.jl
Call ExtractNoiseData!, OrderSlices!, ReverseBipolar!, RemoveRef!.
Check the specific functions for info.
    
MRIReco reference: https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.28792
ISMRMRD reference: https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.26089

# Arguments
* `params::Dict{Symbol, Any}` - MRINavigator parameter structure, check defaultNavParams() for info
"""
function loadRawData(params::Dict{Symbol, Any})

    @info "Load Raw data"
    # load raw data
    rawData = RawAcquisitionData(ISMRMRDFile(params[:path_imgData]),
            repetition = params[:rep])

    if params[:rep] != 0
        for ii = 1:length(rawData.profiles)
            rawData = @set rawData.profiles[ii].head.idx.repetition = 0
        end
    else
        ExtractNoiseData!(rawData) # removing the noise acquisition is only necessary for the first rep
    end
    OrderSlices!(rawData)
    ReverseBipolar!(rawData)
    RemoveRef!(rawData)

    return rawData

end


"""
    convertRawToAcq(rawData::::RawAcquisitionData)

Convert raw data to acquisition data using MRIReco.jl, then apply small adjustments.
Return acquisition data structure.
    
MRIReco reference: https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.28792

# Arguments
* `rawData::RawAcquisitionData` - raw data structure obtained loading raw data with MRIReco.jl
"""
function convertRawToAcq(rawData::RawAcquisitionData)

    @info "convert data and adjust"
    # convert to acquisitionData (note: the estimateProfileCenter flag is set to true)
    acqData = AcquisitionData(rawData, estimateProfileCenter=true)
    CopyTE!(rawData, acqData)
    AdjustSubsampleIndices!(acqData)
    acqData = convertUndersampledData(acqData)

    return acqData

end