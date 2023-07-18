#=
function loadData(params::Dict{Symbol, Any})

    rawData = RawAcquisitionData(ISMRMRDFile(params[:path_imgData]),
                                slice = params[:slices],
                                contrast = params[:echoes],
                                repetition = params[:rep])
    if params[:rep] != 0
        for ii = 1:length(rawData.profiles)
            rawData = @set rawData.profiles[ii].head.idx.repetition = 0
        end
    end
    
    OrderSlices!(rawData)
        if params[:rep] == 0
            noisemat = ExtractNoiseData!(rawData)
            FileIO.save(params[:path_noise],"noisemat",noisemat)
        else
            noisemat = FileIO.load(params[:path_noise], "noisemat")
        end
        ReverseBipolar!(rawData)
        RemoveRef!(rawData, params[:slices], params[:echoes])

        (nav, nav_time) = ExtractNavigator(rawData, params[:slices])
        nav_time = nav_time .* 2.5 # seconds from beginning of the day

    rawMap = RawAcquisitionData(ISMRMRDFile(params[:path_mapData]),
                                slice = params[:slices],
                                contrtast = params[:echoes],
                                repetition = params[:rep])
    


end
=#