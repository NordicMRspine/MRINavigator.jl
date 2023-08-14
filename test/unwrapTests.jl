function test_find_wrapped()

    rawData = FileIO.load(joinpath(datadir, "data.jld2"), "data")
    deleteat!(rawData.profiles, 1:2)
    (nav, nav_time) = ExtractNavigator(rawData)
    nav_time = nav_time[:,1:1] .* 2.5

    # Simulate nav data
    nav_unwrapped = ones(Float64, 1,1,128,1)
    nav_unwrapped[1,1,:,1] = 3 * sin.(Array(0.5:0.5:64)) .+0.2
    nav = angle.(exp.(im * nav_unwrapped))
    wrapped_groundTruth = findall(x -> x > 0.1, abs.(nav .- nav_unwrapped)[1,1,:,:])
    
    # Simulate resp recording
    trace_data = sin.(Array(-5:0.1:68))
    trace_time = range(findmin(nav_time)[1] - 9 * 500, findmax(nav_time)[1] + 8 * 500, length(trace_data))
    trace = hcat(trace_time, trace_data)
    slices = 1
    TR = convert(Int64, rawData.params["TR"]) .* 1e-3

    (wrapped_points, correlation) = find_wrapped(nav, nav_time, trace, slices, TR)
    wrapped_computed = findall(x-> x==1, wrapped_points)

    nav = wrap_corr!(nav, wrapped_points, correlation, 1)

    @test 0.5 < correlation[1] < 0.6
    @test wrapped_computed == wrapped_groundTruth
    
end


function testunwrap()
    @testset "UnwrapTests" begin
        test_find_wrapped()
    end
end

testunwrap()