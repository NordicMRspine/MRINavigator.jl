export ReconstructSaveMap


# FUNCTION TO RECONSTRUCT AND SAVE THE COIL SENSITIVITY MAPS
function ReconstructSaveMap(path_nifti::String, path_ref::String)

    (img, acq) = ReconstructMap(path_ref)
    niftiSaveMap(img, path_nifti, acq)

end


"""
        (img, acq) = econstruct_maps(path_ref::String)

Reconstructs the sensitivity maps.
Returns the image data and acquisition data.
"""
function ReconstructMap(path_ref::String)

    raw = RawAcquisitionData(ISMRMRDFile(path_ref))
    OrderSlice!(raw)
    acq = AcquisitionData(raw, estimateProfileCenter=true)
    sensit = CompSensit(acq)
    img = Reconstruct(acq, sensit)
    
    return img, acq

end

#FUNCTION TO SAVE THE SENSITIVITY MAPS AS NIFTI FILES, FOR COMPATIBILITY WITH SCT
function niftiSaveMap(img::AbstractArray{T,6}, acq::AcquisitionData, path_nifti::String) where {T}

    img = img.data[:, :, :]
    voxel_tmp = fieldOfView(acq)[1:2]./encodingSize(acq)
    voxel_size = (voxel_tmp[1], voxel_tmp[2],fieldOfView(acq)[3])

    orientation = zeros(Float32, 3, 4)
    orientation[1,1] = -voxel_size[1]
    orientation[2,2] = voxel_size[2]
    orientation[3,3] = voxel_size[3]

    ni = NIVolume(abs.(img); voxel_size=voxel_size, orientation = orientation)
    ni = @set ni.header.pixdim[1] = -1.0
    niwrite(path_nifti, ni)

end

