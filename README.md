# MRINavigator
[![Build Status](https://github.com/NordicMRspine/MRINavigator/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/NordicMRspine/MRINavigator/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/NordicMRspine/MRINavigator.jl/graph/badge.svg?token=LMOFFRQIA2)](https://codecov.io/gh/NordicMRspine/MRINavigator.jl)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://NordicMRspine.github.io/MRINavigator.jl/dev)


MRINavigator.jl provides multiple navigator-based correction pipelines for magnetic resonance data. These aim at demodulating time dependent field variations. The package was developed with a focus on spinal cord imaging, however it can be used for multiple imaging applications. 
The corrections are to be applied on the raw data before the image reconstruction. MRIReco.jl can be used to reconstruct the images.

More details can be found in the [online documentation](https://NordicMRspine.github.io/MRINavigator.jl/dev).

Example [data](https://doi.org/10.5281/zenodo.10731729), [scripts and results](https://github.com/NordicMRspine/UserExample_MRINavigator/tree/main) are also available.

MRIReco reference: https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.28792

## How to give credit
If you use this package please acknowledge us by citing:

Beghini L, Büeler S, Liechti MD, Jaffray A, David G, Vannesjo SsJ. [Optimised navigator correction of physiological field fluctuations in multi-echo GRE of the lumbar spinal cord at 3T](https://doi.org/10.1101/2024.12.05.24318389). medRxiv 2024.12.05.24318389; doi: https://doi.org/10.1101/2024.12.05.24318389

## Community Standards
This project is part of the Julia community and follows the [Julia community standards](https://julialang.org/community/standards/).

## Acknowledgement
The authors thank the ISMRM Reproducible Research Study Group for conducting a code review of the code (package version 0.1.1). The scope of the code review covered only the code’s ease of download, quality of documentation, and ability to run, but did not consider scientific accuracy or code efficiency.