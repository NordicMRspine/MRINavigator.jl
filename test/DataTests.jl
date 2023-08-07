
function test_AdjustData_raw(datadir::String)

    rawData = FileIO.load(joinpath(datadir, "data.jld2"), "data")
    flags = ExtractFlags(rawData)

    # test order slices
    number = length(rawData.profiles)
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

    # test reverse bipolar
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

function test_AdjustData_acq(datadir::String)

    rawData = FileIO.load(joinpath(datadir, "data.jld2"), "data")
    acqData = AcquisitionData(rawData, estimateProfileCenter=true)
    CopyTE!(rawData, acqData)
    (nav, nav_time) = ExtractNavigator(rawData)
    selectEcho!(acqData, 0)
    selectSlice!(acqData, 0, nav, nav_time)

    @test acqData.traj[1].TE === convert(typeof(acqData.traj[1].TE), rawData.params["TE"][1])

end


function test_SpineCenterline(datadir::String, tmpResdir::String)
    map = FileIO.load(joinpath(datadir, "map.jld2"), "map")
    OrderSlices!(map)
    acq = AcquisitionData(map, estimateProfileCenter=true)
    sensit = CompSensit(acq)
    img = Reconstruct(acq, sensit)
    niftiSaveMap(img, acq, tmpResdir * "/gre2D_Ref.nii")
    
    @test isfile(tmpResdir * "/gre2D_Ref.nii")

end

function test_CoilSensMap(datadir::String, tmpResdir::String)
    map = FileIO.load(joinpath(datadir, "map.jld2"), "map")
    acqMap = AcquisitionData(map, estimateProfileCenter=true)
    data = FileIO.load(joinpath(datadir, "data.jld2"), "data")
    acqData = AcquisitionData(data, estimateProfileCenter=true)
    sensit = CompSensit(acqMap)
    sensit = ResizeSensit!(sensit, acqMap, acqData)
    img = directreco(acqData)
    sensit_basic = estimateCoilSensitivities(img)
    sensit_basic = sensit_basic[:,:,1,:,:,1]


    err = norm(vec(sensit)-vec(sensit_basic))/norm(vec(sensit_basic))
    @test err < 3

end


function testdata(datadir::String, tmpResdir::String)
    @testset "DataTests" begin
        test_AdjustData_raw()
        test_SpineCenterline()
        test_CoilSensMap()
    end
end