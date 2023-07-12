export defaultNavParams


"""
    params = defaultNavParams()

Define defaul parameters for data loading, navigator correction and image reconstruction.

# Additional required parameters are


# Additional not required parameters are
"""
function defaultNavParams()
    params = Dict{Symbol,Any}()
    params[:slices] = nothing
    params[:echoes] = nothing
    params[:rep] = 0
    params[:reconstruct_map] = true
    params[:comp_sensit] = true
    params[:comp_SCT] = true
    params[:trust_SCT] = false
    params[:corr_type] = "FFT"
    params[:FFT_interval] = 35 # millimiters

  
    return params
  end