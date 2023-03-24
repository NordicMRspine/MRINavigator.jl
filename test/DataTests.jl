
function test_AdjustData(datadir::String)

    map = FileIO.load(joinpath(datadir, "map.jld2"), "map")
    data = FileIO.load(joinpath(datadir, "data.jld2"), "data")
    flags = ExtractFlags(data)
    noisemat = ExtractNoiseData!(data, flags)
    ReverseBipolar!(data, flags)
    flags_Bireverse = ExtractFlags(data)
    params = defaultNavParams()
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
    niftiSavemap(img, acq, tmpResdir)


end


function testdata(datadir::String, tmpResdir::String)
    @testset "DataTests" begin
        test_AdjustData()
        test_SpineCenterline()
    end
end