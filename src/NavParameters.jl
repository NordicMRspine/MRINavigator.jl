export defaultNavParams


"""
    params = defaultNavParams()

Define defaul parameters for data loading, navigator correction and image reconstruction.
The user must add the data paths.
"""
function defaultNavParams()
    params = Dict{Symbol,Any}()
    params[:slices] = nothing
    params[:echoes] = nothing
    params[:rep] = 0
    params[:reconstruct_maps] = true
    params[:comp_sensit] = true
    params[:corr_type] = "FFT"
    params[:FFT_interval] = 35 # millimiters
  
    return params
  end