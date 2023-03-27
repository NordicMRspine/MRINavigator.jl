#=
function loadData(params::Dict{Symbol, Any})

    rawData = RawAcquisitionData(ISMRMRDFile(params[:path_imgData]),
                                repetition = params[:rep])

    raw = RawAcquisitionData(ISMRMRDFile(params[:path_mapData]),
                                slice = params[:slices],
                                contrtast = params[:echoes],
                                repetition = params[:rep])
    


end
=#