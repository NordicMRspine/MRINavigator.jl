# MRINavigator
*Magnetic Resonance Imaging Navigator-based corrections*

## Table of contents

```@contents
Pages = [
    "index.md",
    "GettingStarted.md",
    "Pipelines.md",
    "API.md"
]
Depth = 2
```

## Introduction
MRINavigator provides multiple navigator-based correction pipelines for Magnetic Resonance (MR) images. These aim at demodulating time-dependent field variations present in multi echo-gradient echo acquisitions. The package was developed with a focus on spinal cord imaging, but it can be used for multiple imaging applications. The corrections are to be applied to the raw data before the image reconstruction. [MRIReco.jl](https://github.com/MagneticResonanceImaging/MRIReco.jl) can be used to reconstruct the images.

!!! note
    MRINavigator.jl is newly published, and any feedback is welcome. Please report any bugs or feature requests as an Issue in Github.

## Installation
Start `julia` and open the package manager REPL mode by entering `]`. Then enter
```
pkg> add MRINavigator
```
This will install `MRINavigator` and all its dependencies. If you want to develop
`MRINavigator` itself you can checkout `MRINavigator` locally as usual by calling
```
pkg> dev MRINavigator
```
More information on how to develop a package can be found in the Julia documentation.

### Requirements
To use some package functionalities, external softwares are necessary. 
These include:
* [Spinal Cord Toolbox (SCT)](https://spinalcordtoolbox.com)
* [FSLeyes](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLeyes)
Using these toolboxes should improve the correction outcome of the pipelines including a Fourier transform (FFT) step. Use of SCT and FSLEyes is only relevant for spinal cord acquisitions.
For additional information read the [Get started](@ref) and [Navigator-based correction pipelines](@ref) sections.

## Testing MRINavigator
To make sure that the package is correctly installed, start `julia` from the command line, type `]` to enter the package manager REPL mode. Then enter
```
pkg> test MRINavigator
```

## Updating MRINavigator
To update MRINavigator to the latest version, start `julia` from the command line, type `]` to enter the package manager REPL mode. Then enter
```
pkg> update MRINavigator
```

## Navigator-based correction
Multi-echo gradient-echo (GRE) sequences are commonly acquired both in research and clinical practice. However, one of their main limitations is the sensitivity to field instabilities both in space and time. Indeed, for the signal spatial encoding to be effective, a background homogeneous field in time and space is required. Time-varying background fields can lead to phase modulation between k-space lines, and therefore TE-dependent ghosting artefacts. [Navigator](https://www.sciencedirect.com/science/article/pii/S1053811910003356?via%3Dihub) readouts in the k-space center can be used to measure the intensity of the field fluctuations,enabling correct demodulation of the acquired signal before image reconstruction. The standard navigator-based correction was developed for brain imaging and it is not robust when applied in other areas e.g. the spinal cord. When failing, the correction can even exacerbate the problem. This package provides optimized post-processing pipelines to correct for dynamic field instabilities in GRE sequences. For additional information read the [Navigator-based correction pipelines](@ref) section.

## Plotting
`MRINavigator` does not depend upon a particular plotting package since there are various plotting packages available in Julia. Feel free to use your package of choice. 

## Acknowledgements
This package uses the reconstruction functions and data structures available in [MRIReco.jl](https://github.com/MagneticResonanceImaging/MRIReco.jl).
T. Knopp and M. Grosser (2021). [MRIReco.jl: An MRI Reconstruction Framework written in Julia]( https://doi.org/10.1002/mrm.28792). Magnetic Resonance in Medicine. 2021.

## Citing this work
If you use MRINavigator in you research please cite the following:

[Optimised navigator correction of physiological field fluctuations in multi-echo GRE of the lumbar spinal cord at 3T](https://submissions.mirasmart.com/ISMRM2023/Itinerary/PresentationDetail.aspx?evdid=1673). L Beghini, G David, M D Liechti, S Büeler, S J Vannesjo. 2023. Proceedings of the International Society for Magnetic Resonance in Medicine (ISMRM).