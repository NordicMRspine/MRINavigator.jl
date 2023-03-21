#using MRINavigator
using FileIO
using MRIReco
using Test

using Scratch
using LazyArtifacts

const datadir = joinpath(artifact"TestData")
@info "The test data is located at $datadir."

const tmpdir  = @get_scratch!("tmp")
@info "If you want to check the output of the tests, please head to $tmpdir."


include("DataTests.jl")
