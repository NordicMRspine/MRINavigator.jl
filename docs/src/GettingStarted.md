# Get started

There are available [user examples](https://github.com/NordicMRspine/UserExample_MRINavigator) to get started. The user needs to dispone of their own data to run the pipelines as currently there are not example data available.

## Data acquisition
The user should have the raw data of a gradient echo acquisition. The pipelines can be run on both multi-echo and single-echo data. Other than the main acquisition, which is usually undersampled, also a lower resolution, fully sampled scan is necessary to compute the [coils sennsitivity maps](https://doi.org/10.1002/mrm.1241) and reconstruct the images. This low resolution scan is also called reference scan and should include only one echo. The main gradient echo acquisition must include a navigator readout trough the center of k-space at __the end__ of every TR. During the acqusition it is advisable to connect a repiratory belt and record the signal. This can be used to unwrap the navigator phase estimates if phase wrapping is present.

## Data reshaping
All the data should be exported from the scanner in raw format. Then they should be converted to [ISMRMRD](https://ismrmrd.readthedocs.io/en/latest/index.html) format. Siemens TWIX data can be converted to ISMRMRD using [siemens_to_ismrmrd](https://github.com/ismrmrd/siemens_to_ismrmrd). After the conversion the data can be loaded in the julia framework. The repiratory belt recording must be synchronised with bla bla


## The parameters dictionary
Before calling the package function the correction pipeline should be chosen and the parameters should be set. For more informations regarding the correction pipelines and 


## User examples
* Compact
* Semi-Compact
* Complete

## Disclaimer
MRINavigator and the functions to ajdust the data after loading the were developed using Siemens data [bla bla]