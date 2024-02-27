export defaultNavParams

"""
    params = defaultNavParams()

Define default parameters for data loading, navigator correction and image reconstruction.

# Default parameters options are
* `slices::Union{Nothing, Vector}`    - a vector containing the number of the slices to be loaded, nothing means all slices
* `echoes::Union{Nothing, Vector}`    - a vector containing the number of the echoes to be loaded, nothing means all echoes
* `rep::Int`                          - repetition to be loaded, the first repetition is 0. It is mandatory to select one
* `comp_sensit::Bool`                 - compute the sensitivity maps using the reference scan
* `comp_centerline::Bool`             - use the Spinal Cord Toolbox (SCT) to find the centerlne position
* `trust_SCT::Bool`                   - trust SCT or display the resutls and wait for user feedback with the julia REPL
* `use_centerline::Bool`              - use the spinal cord centerline information in the navigator-based correction
* `corr_type::String`                 - correction type. Options: "none", "knav", "FFT", "FFT_unwrap"
* `FFT_interval::String`              - interval in mm to be considered for the FFT based approach

# Additional required parameters are
* `path_imgData::String`              - path to the image data file in ISMRMRD format
* `path_refData::String`              - path to the reference data file in ISMRMRD format
* `path_sensit::String`               - path to the file where the sensitivity maps will be saved. The file extension must be .mat
* `path_noise::String`                - path to the file where the noise acquisition will be saved. The file extension must be .jld2
* `path_results::String`              - path to the results folder

# Additional optional parameters are
* `path_niftiMap::String`             - path to the file where the reconstructed reference data will be saved in nifti format. The file extension must be .nii
* `path_centerline::String`           - path to the folder where the Spinal Cord Toolbox (SCT) centerline results will be saved
* `path_physio::String`              - path to the physiological trace recording in .mat format. The variable should be a two columns vector (1:time [ms], 2:trace).
                                        The time should be expressed in seconds from the beginning of the day and contain time points before and after the image acquisiton (at least 4 s).

ISMRMRD reference: https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.26089
SCT reference: https://spinalcordtoolbox.com
"""
function defaultNavParams()

    params = Dict{Symbol,Any}()
    params[:slices] = nothing
    params[:echoes] = nothing
    params[:rep] = 0
    params[:comp_sensit] = true
    params[:comp_centerline] = true
    params[:trust_SCT] = false
    params[:use_centerline] = true
    params[:corr_type] = "FFT"
    params[:FFT_interval] = 35 # [millimiters]
  
    return params
  end