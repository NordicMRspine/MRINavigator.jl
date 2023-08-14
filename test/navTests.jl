function test_centerline_position(datadir::String)

    rawMap = FileIO.load(joinpath(datadir, "map.jld2"), "map")
    acqMap = AcquisitionData(rawMap, estimateProfileCenter=true)
    rawData = FileIO.load(joinpath(datadir, "data.jld2"), "data")
    deleteat!(rawData.profiles, 1:2) # remove reference data
    acqData = AcquisitionData(rawData, estimateProfileCenter=true)
    noise = FileIO.load(joinpath(datadir, "noise.jld2"), "noise")
    centerline = [32.0]
    addData = additionalNavInput(noise, rawData, acqData, acqMap, nothing, nothing, centerline)
    centerline = comp_centerline(addData)

    @test centerline == [128]
end

function test_wrap_corr()

    slices = 1
    nav = ones(Float64, 1,1,128,1)
    nav[1,1,:,1] = sin.(Array(1:80:10240))
    wrapped_points = zeros(Int8, 128, 1)
    wrapped_points[10:13,1] .= 1
    # test with positive correlation
    correlation = [0.5]
    nav_unwrapped = deepcopy(nav)
    nav_unwrapped = wrap_corr!(nav_unwrapped, wrapped_points, correlation, slices)

    @test any(nav_unwrapped[1,1,10:13,1] .== nav[1,1,10:13,1] .+ 2pi)

    # test with negaticve correlation
    correlation = [-0.5]
    nav_unwrapped = deepcopy(nav)
    nav_unwrapped = wrap_corr!(nav_unwrapped, wrapped_points, correlation, slices)

    @test any(nav_unwrapped[1,1,10:13,1] .== nav[1,1,10:13,1] .- 2pi)

end

function test_apply_corr(datadir::String)
    
    # Load data
    rawData = FileIO.load(joinpath(datadir, "data.jld2"), "data")
    deleteat!(rawData.profiles, 1:2) # remove reference data
    acqData = AcquisitionData(rawData, estimateProfileCenter=true)
    CopyTE!(rawData, acqData)
    rawMap = FileIO.load(joinpath(datadir, "map.jld2"), "map")
    acqMap = AcquisitionData(rawMap, estimateProfileCenter=true)
    acqData_nocorr = deepcopy(acqData)
    noise = FileIO.load(joinpath(datadir, "noise.jld2"), "noise")
    sensit = CompSensit(acqMap)
    sensit = ResizeSensit!(sensit, acqMap, acqData)

    # Simulate nav data
    nav = ones(Float64, 1,1,128,1)
    nav[1,1,:,1] = 3 .* sin.(Array(1:1:128))

    nav = TE_corr!(nav, acqData, 1e-5, acqData.traj[1].TE .* 1e-3, 256, 1)
    nav = exp.(im .* nav)
    apply_corr!(nav, acqData, 1,128,256,1)

    # Reconstruct the data
    img_corr = Reconstruct(acqData, sensit, noise)
    img = Reconstruct(acqData_nocorr, sensit, noise)

    # Reverse the correction
    nav = ones(Float64, 1,1,128,1)
    nav[1,1,:,1] = - 3 .* sin.(Array(1:1:128))
    nav = TE_corr!(nav, acqData, 1e-5, acqData.traj[1].TE .* 1e-3, 256, 1)
    nav = exp.(im .* nav)
    apply_corr!(nav, acqData, 1,128,256,1)
    img_corrcorr = Reconstruct(acqData, sensit, noise)

    uniformity = measUniformity(img.data, sensit)
    uniformity_corr = measUniformity(img_corr.data, sensit)
    uniformity_corrcorr = measUniformity(img_corrcorr.data, sensit)

    @test uniformity - uniformity_corrcorr < 1e-6
    @test uniformity - uniformity_corr > 1.3

    err = norm(vec(img_corrcorr.data)-vec(img.data))/norm(vec(img.data))
    @test err < 1e-6

    err = norm(vec(img_corr.data)-vec(img.data))/norm(vec(img.data))
    @test err > 1

end

function test_FFTnav_unwrap(datadir::String)
    
    # Load data
    rawMap = FileIO.load(joinpath(datadir, "map.jld2"), "map")
    acqMap = AcquisitionData(rawMap, estimateProfileCenter=true)
    rawData = FileIO.load(joinpath(datadir, "data.jld2"), "data")
    deleteat!(rawData.profiles, 1:2) # remove reference data
    acqData = AcquisitionData(rawData, estimateProfileCenter=true)
    acqData.traj[1].TE = rawData.profiles[2].head.user_int[8] .* 1e-3
    (nav, nav_time) = ExtractNavigator(rawData)
    nav_time = nav_time[:,1:1] .* 2.5
    acqData_nocorr = deepcopy(acqData)
    noise = FileIO.load(joinpath(datadir, "noise.jld2"), "noise")
    sensit = sensit = CompSensit(acqMap)
    sensit = ResizeSensit!(sensit, acqMap, acqData)
    centerline = [32.0]

    # Simulate nav data
    nav = ones(Complex{Float32}, 1,1,128,1)
    nav[1,1,:,1] = exp.(im * 3.25 * sin.(Array(0.5:0.5:64)) .* sin.(Array(0:0.0235:3)))
    nav = repeat(nav, 256, 32, 1, 1)

    # Simulate resp recording
    trace_data = sin.(Array(-4:0.1:68))
    trace_time = range(findmin(nav_time)[1] - 8 * 500, findmax(nav_time)[1] + 9 * 500, length(trace_data))
    trace = hcat(trace_time, trace_data)
    
    # FFT correction
    addData = additionalNavInput(noise, rawData, acqData, acqMap, nav_time, trace, centerline)
    params = defaultNavParams()
    params[:corr_type] = "FFT_unwrap"
    params[:use_SCT] = true
    output = NavCorr!(nav, acqData, params, addData)

    # Reconstruct the data
    img_corr = Reconstruct(acqData, sensit, noise)
    img = Reconstruct(acqData_nocorr, sensit, noise)

    uniformity = measUniformity(img.data, sensit)
    uniformity_corr = measUniformity(img_corr.data, sensit)

    @test uniformity - uniformity_corr > 1.5

    err = norm(vec(img_corr.data)-vec(img.data))/norm(vec(img.data))
    @test err > 0.7

    @test output.centerline == [128]

    @test 0.6 < output.correlation[1] < 0.8

    @test length(findall(x -> x==1, output.wrapped_points)) >= 1
    
end


function testnav(datadir::String)
    @testset "NavigatorTests" begin
        test_centerline_position(datadir)
        test_wrap_corr()
        test_apply_corr(datadir)
        test_FFTnav_unwrap(datadir)
    end
end

testnav(datadir)