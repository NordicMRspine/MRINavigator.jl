
function test_AdjustData_raw(datadir::String)

    rawData = FileIO.load(joinpath(datadir, "data.jld2"), "data")
    flags = ExtractFlags(rawData)

    # test OrderSlices!
    rawData = @set rawData.profiles[1].head.position[3] = rawData.profiles[1].head.position[3] + 1
    rawData = @set rawData.profiles[1].head.idx.slice = 2
    position = rawData.profiles[1].head.idx.slice
    OrderSlices!(rawData)
    position_ordered = rawData.profiles[1].head.idx.slice

    @test position_ordered < position

    # the noise acquision has flag 19
    rawData.profiles[1].head.flags = rawData.profiles[1].head.flags + 2^18
    noisemat_rawData = rawData.profiles[1].data
    noisemat = ExtractNoiseData!(rawData)

    @test any(flags_Bireverse[:,19] .== false)
    @test noisemat == noisemat_rawData

    # test reverse bipolar. The revese flag is number 22
    rawData.profiles[1].head.flags = rawData.profiles[1].head.flags + 2^21
    reversed_profile = rawData.profiles[1].data
    ReverseBipolar!(rawData)
    flags_Bireverse = ExtractFlags(rawData)

    @test any(flags_Bireverse[:,22] .== false)
    @test reversed_profile == reverse!(rawData.profiles[1].data)
    
    # check number of profiles
    numflags = size(flags,1)
    numProfiles = size(rawData.profiles, 1)
    RemoveRef!(rawData)
    slices = rawData.params["enc_lim_slice"].maximum + 1
    echoes = size(rawData.params["TE"],1) + 1
    
    @test size(rawData.profiles, 1) == numProfiles - (slices * echoes)
    @test numProfiles == numflags - 1

end

function test_AdjustData_acq(datadir::String, tmpResdir::String)

    numProfiles = size(rawData.profiles, 1)
    rawData = FileIO.load(joinpath(datadir, "data.jld2"), "data")
    acqData = AcquisitionData(rawData, estimateProfileCenter=true)
    CopyTE!(rawData, acqData)
    (nav, nav_time) = ExtractNavigator(rawData)
    nav_rawData = zeros(ComplexF32, size(rawData.profiles[1].data)[1], size(rawData.profiles[1].data)[2],
        rawData.params["reconSize"][2], rawData.params["enc_lim_slice"].maximum +1)
    for ii = 4:2:numProfiles
        nav_rawData[:,:,div(ii-2,2),1] = rawData.profiles[ii].data
    end

    @test acqData.traj[1].TE == convert(typeof(acqData.traj[1].TE), rawData.params["TE"][1])
    # The first 2 profiles are reference data
    @test nav[:,:,:,:] == nav_rawData
    @test nav_time[1,1] == convert(typeof(nav_time[1,1]), rawData.profiles[4].head.acquisition_time_stamp)

end


function test_CoilSensMap(datadir::String, tmpResdir::String)

    map = FileIO.load(joinpath(datadir, "map.jld2"), "map")
    acqMap = AcquisitionData(map, estimateProfileCenter=true)
    data = FileIO.load(joinpath(datadir, "data.jld2"), "data")
    acqData = AcquisitionData(data, estimateProfileCenter=true)
    sensit = CompSensit(acqMap)
    #binarize sensit
    thresh = 0.5* mean(abs.(sensit))
    cartes_index_binar = findall(x -> x > thresh, abs.(sensit))
    sensit_binar = zeros(Int64, size(sensit))
    sensit_binar[cartes_index_binar] .= 1
    sensit_binar = circshift(sensit_binar, (0,-1,0,0))
    #resize sensit
    sensit_resized = ResizeSensit!(sensit, acqMap, acqData)
    #compare sensit with basic version
    img = directreco(acqData)
    sensit_basic = estimateCoilSensitivities(img)
    sensit_basic = sensit_basic[:,:,1,:,:,1]

    err = norm(vec(sensit_resized)-vec(sensit_basic))/norm(vec(sensit_basic))
    # test whole algorithm
    @test err < 3
    # test ResizeSensit!
    @test size(sensit_resized) == size(sensit_basic)
    # test removeBehindBack!
    err = norm(vec(reverse(sensit_binar[:,33:end,1,:], dims = 2))-vec(sensit_binar[:,1:32,1,:]))/norm(vec(sensit_binar[:,1:32,1,:]))
    @test err < 0.15

    # Save sensitivity maps in a temporary folder
    FileIO.save(joinpath(tmpResdir, "sensit.jld2"), "sensit", sensit)

end

function test_niftsave(datadir::String, tmpResdir::String)

    map = FileIO.load(joinpath(datadir, "map.jld2"), "map")
    acq = AcquisitionData(map, estimateProfileCenter=true)
    sensit = FileIO.load(joinpath(tmpResdir, "sensit.jld2"), "sensit")
    img = Reconstruct(acq, sensit)
    niftiSaveImg(img, acq, tmpResdir * "/gre2D_Ref.nii")
    
    @test isfile(tmpResdir * "/gre2D_Ref.nii")

end


function testdata(datadir::String, tmpResdir::String)
    @testset "DataTests" begin
        test_AdjustData_raw()
        test_CoilSensMap()
        test_niftisave()
    end
end