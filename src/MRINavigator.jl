module MRINavigator

using MRIReco
using MRICoilSensitivities
using MRIFiles
using Statistics
using Setfield
using Images
using PolygonOps
using NIfTI
using REPL.TerminalMenus
using DataInterpolations
using DSP
using FileIO
using MAT
using CSV
using DataFrames
using FFTW


include("AdjustData.jl")
include("CoilSensMap.jl")
include("main.jl")
include("NavData.jl")
include("NavParameters.jl")
include("Reconstruction.jl")
include("SpineCenterline.jl")
include("Navigator.jl")
include("Unwrap.jl")

end # module
