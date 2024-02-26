# Get started

There are available [user examples](https://github.com/NordicMRspine/UserExample_MRINavigator) to get started.
Example data can be downloaded [here](). The user can also test the pipelines on their own data.

## Data requirements
The navigator based correction can only be applied on the raw data of a gradient echo acquisitions in [MRD format](https://ismrmrd.readthedocs.io/en/latest/index.html). The gradient echo acquisition can be multi-echo or single-echo. One requirement on the acquisiton is to have maximum one concatenation, meaning that the sequence repetition time (TR) shold be long enough to acquire one line in every slice. This is necessary to correctly identify and remove the reference data from the Siemens scans. The gradient echo acquisition must include a navigator readout through the center of k-space at __the end__ of each TR. During the acquisition, it is advisable to collect the signal from a respiratory belt as a reference. This can be used to unwrap the navigator's phase estimates if phase wrapping is present.

Other than the main acquisition, which is usually undersampled, a lower resolution, fully sampled scan is also necessary to compute the [coils sensitivity maps](https://doi.org/10.1002/mrm.1241) and reconstruct the images. This low-resolution scan is also called a reference scan and should include only one echo. It is possible to extract a single echo from a multi echo acquisition using the `selectEcho!` function. 

## Data reshaping
All the data should be exported from the scanner in raw format. Then they should be converted to [ISMRMRD](https://ismrmrd.readthedocs.io/en/latest/index.html) format. Siemens TWIX data can be converted to ISMRMRD using [siemens_to_ismrmrd](https://github.com/ismrmrd/siemens_to_ismrmrd). After the conversion, the data can be loaded into the Julia framework. Conversion of data from other vendors has not been explicitly tested by the authors.
The repiratory belt recording must be synchronised with the time stamps in the image acquisition (i.e resampled). Then they must be saved in a two-column vector (1:time [ms], 2:trace) in .mat format. Each repetition should be in a different file. The time should be expressed in seconds from the beginning of the day and contain time points before and after the image acquisition (at least 4 s).

## The parameters dictionary
Before calling the package functions, the relevant correction pipeline should be chosen and the parameters dictionary should be filled. Also the data paths and results paths need to be defined. For more details regarding the correction pipelines and parameters read the [Navigator-based correction pipelines](@ref) page.
All the information necessary to apply the corrections is defined in a [dictionary](https://docs.julialang.org/en/v1/base/collections/#Dictionaries). This includes all the file paths and analysis parameters. The user can also add items to the dictionary if needed.
Here is an example of a `params` dictionary:

```julia
params = Dict{Symbol,Any}()
params[:subject] = "sub_01"
params[:slices] = [1,2] # type nothing for all slices
params[:echoes] = [3,4] # type nothing for all echoes
params[:rep] = 0
params[:comp_sensit] = true
params[:comp_centerline] = true
params[:trust_SCT] = false
params[:use_centerline] = true
params[:corr_type] = "FFT_unwrap"
params[:FFT_interval] = 35 # millimetres
params[:root_path] = "/Users/me/my_data/"

params[:label] = params[:corr_type] * "_rep_" * string(params[:rep])
params[:path_imgData] = params[:root_path] * params[:subject] * "/h5/gre2D.h5"
params[:path_refData] = params[:root_path] * params[:subject] * "/h5/gre2D_Ref.h5"
params[:path_niftiMap] = params[:root_path] * params[:subject] * "/Nifti/gre2D_Ref.nii"
params[:path_centerline] = params[:root_path] * params[:subject] * "/Nifti/"
params[:path_physio] = params[:root_path] * params[:subject] * "/Physiological_trace/belt_reco_rep"
params[:path_sensit] = params[:root_path] * params[:subject] * "/Results/senseMap_GRE.jld2"
params[:path_noise] = params[:root_path] * params[:subject] * "/Results/noisemat.jld2"
params[:path_results] = params[:root_path] * params[:subject] * "/Results/"
params[:file_name] = "gre2D"
```

## User examples
Three user examples are available in the folder [user examples](https://github.com/NordicMRspine/UserExample_MRINavigator):
* __Compact__: runs all the selected pipeline automatically but it is not customizable and not amenable to debugging.
* __Semi-Compact__: allows for some level of customization and it is easy to debug.
* __Complete__: requires more knowledge of the data structures but it is flexible and adaptable.

## Disclaimer
Siemens data only were used to develop MRINavigator. All the functions to adjust the data before running the pipeline (e.g., the function to extract the navigator profiles or to remove the reference profiles) have been tested on Siemens data only. There is no guarantee that all of these functions are needed and will work on other vendors data. Other vendors users should convert the raw data in [ISMRMRD](https://ismrmrd.readthedocs.io/en/latest/index.html) format and when loading these into the Julia framework they should make sure that all the needed information is present. Please start from the [complete user example](https://github.com/NordicMRspine/UserExample_MRINavigator) if doing this. The functions to compute and apply the corrections should then work correctly.