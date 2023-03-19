module MRINavigator

using MRIReco
using MRICoilSensitivitiess


include("AdjustData.jl")
include("CoilSensMap.jl")
include("NavParameters.jl")
include("Reconstruction.jl")
include("SpineCenterline.jl")



end # module
