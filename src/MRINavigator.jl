module MRINavigator

using MRIReco
using MRICoilSensitivities
using MRIFiles
using Setfield
using Images
using PolygonOps
using NIfTI

include("AdjustData.jl")
include("CoilSensMap.jl")
include("SetupData.jl")
include("NavData.jl")
include("NavParameters.jl")
include("Reconstruction.jl")
include("SpineCenterline.jl")

end # module
