module MRINavigator

using MRIReco
using MRICoilSensitivities
using MRIFiles
using Setfield
using Images
using PolygonOps
using NIfTI
using REPL.TerminalMenus

include("AdjustData.jl")
include("CoilSensMap.jl")
include("main.jl")
include("NavData.jl")
include("NavParameters.jl")
include("Reconstruction.jl")
include("SpineCenterline.jl")
#include("Navigator.jl")
include("Unwrap.jl")

end # module
