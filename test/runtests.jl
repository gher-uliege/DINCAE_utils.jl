using Test
using DINCAE_utils

@testset "linreg" begin

    X  = randn(10000,3); b = [1,2,3]; y = 1 .+ X*b

    a2,b2,R22 = DINCAE_utils.linreg(X,y)
    @test a2 ≈ 1
    @test b2 ≈ b
    @test R22 ≈ 1
end
