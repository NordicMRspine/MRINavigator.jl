
function measUniformity(img::Matrix{T}, sensit::Array{T,4}) where {T}

    maxval = findmax(abs.(img))[1]
    img = img ./ maxval
    thresh = 0.5* mean(abs.(sensit[:,:,1,1]))
    mask_index = findall(x -> x > thresh, abs.(sensit[:,:,1,1]))
    uniformity = 1 / std(abs.(img[mask_index]))

    return uniformity

end