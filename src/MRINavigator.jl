module MRINavigator

using MRIReco
using MRICoilSensitivities
using Setfield
using Images
using PolygonOps

include("AdjustData.jl")
include("CoilSensMap.jl")
include("LoadData.jl")
include("NavData.jl")
include("NavParameters.jl")
include("Reconstruction.jl")
include("SpineCenterline.jl")

end # module
