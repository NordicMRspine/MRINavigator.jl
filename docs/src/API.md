# API

This page contains documentation of the public API of MRINavigator. In the Julia REPL one can access this documentation by entering the help mode with `?` and then writing the function for which the documentation should be shown.

# Run compact pipeline
```@docs
MRINavigator.defaultNavParams
MRINavigator.runNavPipeline
MRINavigator.saveNoise
MRINavigator.loadRawData
MRINavigator.convertRawToAcq
```


## Coil sensitivity maps
```@docs
MRINavigator.CompSensit
MRINavigator.CompRoughMask
MRINavigator.ResizeSensit!
MRINavigator.CompResizeSaveSensit
```

## Find centerline
```@docs
MRINavigator.findCenterline
MRINavigator.ReconstructMap
MRINavigator.ReconstructSaveMap
MRINavigator.callSCT
MRINavigator.comp_centerline_pos
```

## Utils
```@docs
MRINavigator.Reconstruct
MRINavigator.directreco
MRINavigator.niftiSaveImg
```

## Navigator correction
```@docs
MRINavigator.NavCorr!
MRINavigator.wrap_corr!
MRINavigator.find_wrapped
MRINavigator.TE_corr!
MRINavigator.apply_corr!
```

## Adjust data
```@docs
 MRINavigator.OrderSlices!
 MRINavigator.ExtractFlags
 MRINavigator.ExtractNoiseData!
 MRINavigator.ReverseBipolar!
 MRINavigator.RemoveRef!
 MRINavigator.CopyTE!
 MRINavigator.AdjustSubsampleIndices!
 MRINavigator.ExtractNavigator
 MRINavigator.selectEcho!
 MRINavigator.selectSlice!
 MRINavigator.additionalNavInput
```