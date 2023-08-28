# MRINavigator
*Magnetic Resonance Imaging Navigator-based corrections*

## Table of contents

```@contents
Pages = [
    "index.md",
    "Pipelines.md",
    "GettingStarted.md",
    "API.md"
]
Depth = 2
```

## Introduction
MRINavigator provides multiple navigator-based correction pipelines for Magnetic Resonance (MR) images. These aim at demodulating time-dependent field variations present in multi echo-gradient echo acquisitions. The package was developed with a focus on spinal cord imaging, however it can be used for multiple imaging applications. The corrections are to be applied on the raw data before the image reconstruction. [MRIReco.jl](https://github.com/MagneticResonanceImaging/MRIReco.jl) can be used to reconstruct the images.

!!! note
    MRINavigator.jl is newly published and any feedback is more than welcome.

## Installation
Start `julia` and open the package manager REPL mode by entering `]`. Then enter
```julia
pkg> add MRINavigator
```
This will install `MRINavigator` and all its dependencies. If you want to develop
`MRINavigator` itself you can checkout `MRINavigator` by calling
```julia
pkg> dev MRINavigator
```
More information on how to develop a package can be found in the Julia documentation.

### requirements
To use some package functionalities external softwares are necessary. These include
* [Spinal Cord Toolbox (SCT)](https://spinalcordtoolbox.com)
* [FSLeyes](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLeyes)
Using these should improve the correction outcome of the pipelines including a fourier transform (FFT) step. This is only relevant for spinal cord acquisitions.
For additinal information read the [Getting Started](@ref) and [Pipelines](@ref) sections.

## Testing MRINavigator
To make sure that the package is correctly installed and works, start `julia` from the command line, type `]` to enter the package manager REPL mode. Then enter
```julia
pkg> test MRINavigator
```

## Updating MRINavigator
To update MRINavigator to the latest version, start `julia` from the command line, type `]` to enter the package manager REPL mode. Then enter
```julia
pkg> update MRINavigator
```
## Navigator-based correction
Multi-echo gradient-echo (GRE) sequence are commonly acquired both in research labs and clinical practice. However, one of their main limitations is the sensitivity to field instabilities both in space and time. Indeed, for the signal spatial encoding to be effective a background homogeneous field in time and space is required. Time-varying background fields can lead to phase modulation between k-space lines, and therefore TE-dependent ghosting artefacts. [Navigator](https://www.sciencedirect.com/science/article/pii/S1053811910003356?via%3Dihub) readouts in the k-space center can be used to measure the intensity of the filed fluctuations allowing to demodulate the acquired signal before the image reconstruction. The standard navigaotr-based correction was developed for brain imaging and it is not robus when applied in other areas e.g. the spinal cord. When failing the correction can even exhacerbate the artifacts. This package provides optimised post-processing pipelines to correct for dynamic field instabilities in GRE sequences. For additinal information read the [Pipelines](@ref) section.

## Plotting
`MRINavigator` is not depending on a particular plotting package since there
are various plotting packages available in Julia.

## Citing this work
If you use MRINavigator in you research please cite the following:

[Optimised navigator correction of physiological field fluctuations in multi-echo GRE of the lumbar spinal cord at 3T](https://submissions.mirasmart.com/ISMRM2023/Itinerary/PresentationDetail.aspx?evdid=1673). L Beghini, G David, M D Liechti, S Büeler, S J Vannesjo. 2023. Proceedings of the International Society for Magnetic Resonance in Medicine (ISMRM).