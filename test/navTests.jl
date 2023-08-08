
function test_knav(datadir::String, tmpResdir::String)
    
    # Load data
    rawMap = FileIO.load(joinpath(datadir, "map.jld2"), "map")
    acqMap = AcquisitionData(map, estimateProfileCenter=true)
    rawData = FileIO.load(joinpath(datadir, "data.jld2"), "data")
    deleteat!(rawData.profiles, 1:2) # remove reference data
    acqData = AcquisitionData(data, estimateProfileCenter=true)
    noise = FileIO.load(joinpath(datadir, "noise.jld2"), "noise")
    sensit = FileIO.load(joinpath(tmpResdir, "sensit.jld2"), "sensit")
    sensit = ResizeSensit!(sensit, acqMap, acqData)

    # Simulate nav data

    # Simulate centerline

    # Simulate trac

    # Navigator correction

    # Reconstruct the data
    
    img = Reconstruct(acqData, sensit, noisemat)

    @test

end


function test(datadir::String, tmpResdir::String)
    @testset "NavigatorTests" begin
        test_knav()
        test_navunwrap()
    end
end