export defaultNavParams

"""
    params = defaultNavParams()

Define defaul parameters for data loading, navigator correction and image reconstruction.

# Additional required parameters are
* `path_imgData::String`              - path to the image data file in ISMRMRD format
* `path_refData::String`              - path to the reference data file in ISMRMRD format
* `path_sensit::String`               - path to the file where the sensitivity maps will be saved. The file extension must be .mat
* `path_noise::String`                - path to the file where the noise acquisition will be saved. The file extension must be .jld2
* `path_results::String`              - path to the results folder

# Additional not required parameters are
* `path_niftiMap::String`             - path to the file where the reconstructed reference data will be saved in nifti format. The file extension must be .nii
* `path_centerline::String`           - path to the folder where the Spinal Cord Toolbox (SCT) centerline results will be saved
* `params[:path_physio]`              - path to the physiological trace recording in .mat format

ISMRMRD reference: https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.26089
SCT reference: https://spinalcordtoolbox.com

"""
function defaultNavParams()
    params = Dict{Symbol,Any}()
    params[:slices] = nothing
    params[:echoes] = nothing
    params[:rep] = 0
    params[:reconstruct_map] = false
    params[:comp_sensit] = true
    params[:comp_SCT] = false
    params[:trust_SCT] = false
    params[:corr_type] = "FFT"
    params[:FFT_interval] = 35 # millimiters

  
    return params
  end