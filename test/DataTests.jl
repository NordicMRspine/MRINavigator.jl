
function test_niftisavemap(path = "test/Data/")

    map = FileIO.load(path, "map")

end


function testdata(path = "test/Data/")
    @testset "DataTests" begin
        test_niftisavemap()
    end
end