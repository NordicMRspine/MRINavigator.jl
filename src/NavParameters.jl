export defaultNavParams

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