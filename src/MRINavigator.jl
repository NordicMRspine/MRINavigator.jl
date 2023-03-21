module MRINavigator

using MRIReco
using MRICoilSensitivities
using Setfield


include("AdjustData.jl")
include("CoilSensMap.jl")
include("NavParameters.jl")
include("Reconstruction.jl")
include("SpineCenterline.jl")



end # module
