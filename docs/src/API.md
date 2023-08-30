# API

This page contains documentation of the public API of MRINavigator. In the Julia REPL one can access this documentation by entering the help mode with ? and then writing the function for which the documentation should be shown. 
For example: `? findCenterline`

# Run compact pipeline
```@docs
MRINavigator.defaultNavParams :: Tuple{}
MRINavigator.runNavPipeline :: Dict{Symbol, Any}
MRINavigator.saveNoise :: Tuple{String, String}
MRINavigator.loadRawData :: Tuple{Dict{Symbol, Any}}
MRINavigator.convertRawToAcq :: Tuple{MRIBase.RawAcquisitionData}
```

<!--
## Coil sensitivity maps
```@docs
MRINavigator.CompSensit :: Union{Tuple{MRIBase.AcquisitionData}, Tuple{MRIBase.AcquisitionData, Any}}
MRINavigator.CompRoughMask :: Tuple{MRIBase.AcquisitionData, Int64, Any}
MRINavigator.ResizeSensit! :: Union{Tuple{T}, Tuple{Array{Complex{T}, 4}, MRIBase.AcquisitionData, MRIBase.AcquisitionData}} where T
MRINavigator.CompResizeSaveSensit :: Union{Tuple{MRIBase.AcquisitionData, MRIBase.AcquisitionData, String}, Tuple{MRIBase.AcquisitionData, MRIBase.AcquisitionData, String, Any}}
```

## Find centerline
```@docs
MRINavigator.findCenterline :: Tuple{Dict{Symbol, Any}}
MRINavigator.ReconstructMap :: Tuple{String}
MRINavigator.ReconstructSaveMap :: Tuple{String, String}
MRINavigator.callSCT :: Tuple{Dict{Symbol, Any}}
MRINavigator.comp_centerline_pos :: Tuple{additionalNavInput}
```

## Utils
```@docs
MRINavigator.Reconstruct :: AcquisitionData, Array{Complex{T},4}, Union{Array{Complex{T}},Nothing} where {T}
MRINavigator.directreco :: Tuple{MRIBase.AcquisitionData}
MRINavigator.niftiSaveImg :: Union{Tuple{T}, Tuple{AbstractArray{T}, MRIBase.AcquisitionData, String}} where T
```

## Navigator correction
```@docs
MRINavigator.NavCorr! :: Array{Complex{T}, 4}, AcquisitionData, Dict{Symbol, Any}, additionalNavInput where{T}
MRINavigator.wrap_corr! :: Tuple{Array{Float64, 4}, Matrix{Int8}, VecOrMat{Float64}, Int64}
MRINavigator.find_wrapped :: Tuple{Array{Float64, 4}, Matrix{Float64}, Matrix{Float64}, Int64}
MRINavigator.TE_corr! :: Union{Tuple{T}, Tuple{Array{T, 4}, MRIBase.AcquisitionData, Float64, Float64, Int64, Int64}} where T
MRINavigator.apply_corr! :: Union{Tuple{T}, Tuple{Array{T, 4}, MRIBase.AcquisitionData, Int64, Int64, Int64, Int64}} where T
```

## Adjust data
```@docs
 MRINavigator.OrderSlices! :: Tuple{MRIBase.RawAcquisitionData}
 MRINavigator.ExtractFlags :: Tuple{MRIBase.RawAcquisitionData}
 MRINavigator.ExtractNoiseData! :: Tuple{MRIBase.RawAcquisitionData}
 MRINavigator.ReverseBipolar! :: Tuple{MRIBase.RawAcquisitionData}
 MRINavigator.RemoveRef! :: Tuple{MRIBase.RawAcquisitionData}
 MRINavigator.CopyTE! :: Tuple{MRIBase.RawAcquisitionData, MRIBase.AcquisitionData}
 MRINavigator.AdjustSubsampleIndices! :: Tuple{MRIBase.AcquisitionData}
 MRINavigator.ExtractNavigator :: Tuple{MRIBase.RawAcquisitionData}
 MRINavigator.selectEcho! :: Tuple{MRIBase.AcquisitionData, Vector{Int64}}
 MRINavigator.selectSlice! :: Union{Tuple{T}, Tuple{MRIBase.AcquisitionData, Vector{Int64}}, Tuple{MRIBase.AcquisitionData, Vector{Int64}, Union{Nothing, Array{Complex{T}, 4}}}, Tuple{MRIBase.AcquisitionData, Vector{Int64}, Union{Nothing, Array{Complex{T}, 4}}, Union{Nothing, Matrix{Float64}}}} where T
 MRINavigator.additionalNavInput :: Union{Tuple{Matrix{ComplexF32}, MRIBase.RawAcquisitionData, MRIBase.AcquisitionData}, Tuple{Matrix{ComplexF32}, MRIBase.RawAcquisitionData, MRIBase.AcquisitionData, Union{Nothing, MRIBase.AcquisitionData}}, Tuple{Matrix{ComplexF32}, MRIBase.RawAcquisitionData, MRIBase.AcquisitionData, Union{Nothing, MRIBase.AcquisitionData}, Union{Nothing, Matrix{Float64}}}, Tuple{Matrix{ComplexF32}, MRIBase.RawAcquisitionData, MRIBase.AcquisitionData, Union{Nothing, MRIBase.AcquisitionData}, Union{Nothing, Matrix{Float64}}, Union{Nothing, Matrix{Float64}}}, Tuple{Matrix{ComplexF32}, MRIBase.RawAcquisitionData, MRIBase.AcquisitionData, Union{Nothing, MRIBase.AcquisitionData}, Union{Nothing, Matrix{Float64}}, Union{Nothing, Matrix{Float64}}, Union{Nothing, Vector{Float64}}}}
```

-->