
function test_AdjustData(datadir::String)

    data = FileIO.load(joinpath(datadir, "data.jld2"), "data")
    flags = ExtractFlags(data)
    noisemat = ExtractNoiseData!(data, flags)
    ReverseBipolar!(data, flags)
    flags_Bireverse = ExtractFlags(data)
    numProfiles = size(data.profiles, 1)
    RemoveRef!(data, 1, 1)

    @test any(flags_Bireverse[:,22] .== false)
    @test any(flags_Bireverse[:,19] .== false)
    @test size(data.profiles, 1) == numProfiles - 2

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
        test_AdjustData()
        test_SpineCenterline()
        test_CoilSensMap()
    end
end