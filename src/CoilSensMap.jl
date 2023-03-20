# FUNCTION TO COMPUTE THE COIL SENSITIVITY MAP
function CompSensit(acq::AcquisitionData)

    sensit = espirit(acq,(6,6),30,eigThresh_1=0.005, eigThresh_2=0)
    slices = numSlices(acq)
    coils = numChannels(acq)
    # compute mask
    mask = CompRoughMask(acq, slices)

    for ii = 1:slices
        components = label_components(mask[:,:,ii])
        measured_area = component_lengths(components)
        measured_area = measured_area[2:end] #remove background component
        blob = findmax(measured_area)[2]
        cartes_index_blob = findall(x -> x!=blob, components)
        mask_slice[cartes_index_blob, ii] .= 0

        # remove noisy voxels on the left of the image
        # corresponding to the back of the subject
        density = sum(mask_slice, dims= 1)[1,:]
        # compute dervative of points density
        lines = div(size(mask_slice, 2),2)
        dder = zeros(Int64, lines)
        for ii = 1:lines
            dder[ii] = density[ii+1] - density[ii]
        end
        # put to zero everything behind the back
        position = findmax(dder)[2] -3
        for jj=1:position
            mask_slice[:,jj].=0
        end
        
        #the convex_hull function requires a vector of vetors in imput containing the points coordinates
        convex_hull_points_old = Tuple.(findall(x -> x==1, mask_slice))
        points_num = size(convex_hull_points_old)[1]
        convex_hull_points = [Vector{Int64}(undef, 2) for _=1:points_num]
        for ii = 1:points_num
            convex_hull_points[ii] = collect(convex_hull_points_old[ii])
        end
        hull = convex_hull(convex_hull_points)
        push!(hull, hull[1])
        cartes_index_slice_old1 = findall(x -> x>(-1), s[:,:,ii])
        cartes_index_slice_old = Tuple.(cartes_index_slice_old1)
        points_num = size(cartes_index_slice_old)[1]
        cartes_index_slice = [Vector{Int64}(undef, 2) for _=1:points_num]
        for ii = 1:points_num
            cartes_index_slice[ii] = collect(cartes_index_slice_old[ii])
        end
        inside = [inpolygon(p, hull; in=true, on=true, out=false) for p in cartes_index_slice]
        #plot!(p, VPolygon(hull), alpha=0.2)
        inside = .!(inside)
        deleteat!(cartes_index_slice_old1, inside);
        mask_slice = mask[:,:,ii]
        mask_slice .= 0
        mask_slice[cartes_index_slice_old1] .= 1
        mask[:,:,ii] = mask_slice
    end
    for ii=1:coils
        sensit[:,:,:,ii] = sensit[:,:,:,ii] .* mask
    end

    return sensit

end

function CompRoughMask(acq::AcquisitionData, slices::Int64)

    img = directreco(acq)
    thresh = 0.14
    I_sum = sqrt.(sum(abs.(img) .^ 2, dims = 5)) .+ eps()
    I_sum = dropdims(I_sum, dims = tuple(findall(size(I_sum) .== 1)...))
    I_max = ones(Float64, slices)
    mask = zeros(size(I_sum))
    for ii = 1:slices
        I_max[ii] = maximum(abs.(I_sum[:,:,ii]))
        mask[findall(x -> x > thresh * I_max[ii], I_sum[:,:,ii]), ii] .= 1
    end

    return mask
end