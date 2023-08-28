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
MRINavigator provides tools for applying a navigator-based correction to Magnetic Resonance Images.

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

## Updating MRINavigator
To update MRINavigator to the latest version, start `julia` from the command line, type `]` to enter the package manager REPL mode. Then enter
```julia
pkg> update MRINavigator
```
## Navigator-based correction

## Plotting
`MRINavigator` is not depending on a particular plotting package since there
are various plotting packages available in Julia.

## Citing this work
If you use MRINavigator in you research please cite the following:

To be published