using PGA
using Test

@testset "PGA.jl" begin
    # Write your tests here.
    @test 1 + 1 == 2
    @test see(2,3) == 5
end
