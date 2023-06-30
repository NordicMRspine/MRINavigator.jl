export ReconstructSaveMap, ReconstructMap, niftiSaveImg, callSCT


"""
    ReconstructSaveMap(path_nifti::String, path_ref::String)

Reconstruct the coil sensitivity map and save it in nifti format without spatial informations.

# Arguments
* `path_nifti::String` - path of the nifti file
* `path_rep::String` - path of reference data in ISMRMRD format
"""
function ReconstructSaveMap(path_nifti::String, path_ref::String)

    (img, acq) = ReconstructMap(path_ref)
    start_voxel = div(acq.encodingSize[1] -  acq.encodingSize[2], 2)
    img = img[start_voxel+1 : start_voxel+acq.encodingSize[2], :, :]
    niftiSaveImg(img, acq, path_nifti)

end

function ReconstructMap(path_ref::String)

    raw = RawAcquisitionData(ISMRMRDFile(path_ref))
    OrderSlices!(raw)
    acq = AcquisitionData(raw, estimateProfileCenter=true)
    sensit = CompSensit(acq)
    img = Reconstruct(acq, sensit)
    
    return img, acq

end


"""
    niftiSaveImg(img::AbstractArray{T}, acq::AcquisitionData, path_nifti::String)

Save the module of the reconstruction output in nifti format, without spatial information.

# Arguments
* `img::AbstractArray{T}` - reconstruction output
* `acq::AcquisitionData` - reconstruction input needed for saving the voxel dimension
* `path_nifti::String` - path of the nifti file
"""
function niftiSaveImg(img::AbstractArray{T}, acq::AcquisitionData, path_nifti::String) where {T}

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


"""
    callSCT(params::Dict{Symbol, Any})

Call spinal cord toolbox and find spinal cord centerline.
https://spinalcordtoolbox.com

# Arguments
* `params::Dict{Symbol, Any}` - paramerters dictionary
"""
function callSCT(params::Dict{Symbol, Any})

    path_nifti = params[:path_niftiMap]
    path_centerline = params[:path_centerline] * "centerline.nii"
    run(`sct_get_centerline -i $path_nifti -c t2s -o $path_centerline`)
    if params[:trust_SCT] == false
        run(`fsleyes $path_nifti -cm greyscale $path_centerline -cm red`)
        options = ["yes", "no"]
        menu = RadioMenu(options)
        choice = request("Is the result satisfying?", menu)
        if choice == 2
            run(`sct_get_centerline -i $path_nifti -c t2s -o $path_centerline -method viewer`)
        end
    end
end