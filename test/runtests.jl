using MRINavigator
using FileIO
using Test
using Coverage
using Setfield

using Scratch
using LazyArtifacts

const datadir = joinpath(artifact"TestDataNavigator", "data")
@info "The test data is located at $datadir."

const tmpResdir  = @get_scratch!("tmp")
@info "If you want to check the output of the tests, please head to $tmpResdir."

include("testFunction.jl")
include("dataTests.jl")
include("navTests.jl")
include("unwrapTests.jl")