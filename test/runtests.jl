using Revise
using Test

using Base.Threads; println("--> There are ", Threads.nthreads(), " threads")

@testset "BifurcationKit" begin
    include("bratu.jl")
end