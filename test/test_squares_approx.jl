# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).

using RunStatistics
using Test

# I measured the timing for N=96 with the exact and approximate solution. See if julia is really faster. It seems some experiment with threading was done but without threads, I'd be surprised if it's faster than the C++ code built in release mode

@testset "squares_approx" begin

    T_obs = 20
    L = 800
    N = 80
    n = 10
    Ns = [N, n]
    epsp_1 = 10^(-2)
    epsp_2 = 0
        

    @test RunStatistics.squares_pvalue_approx(T_obs, Ns, epsp_1) ≈ 0.07082230261169509
    @test RunStatistics.squares_cdf_approx(T_obs, Ns, epsp_1) ≈ 0.9291776973883049

    @test RunStatistics.squares_cdf_approx(T_obs, L, epsp_1) ≈ 0.9291776973883049
    @test RunStatistics.squares_pvalue_approx(T_obs, L, epsp_1) ≈ 0.07082230261169509


    @test RunStatistics.squares_pvalue_approx(T_obs, Ns, epsp_2) == 0.0708223271587789
    @test RunStatistics.squares_cdf_approx(T_obs, Ns, epsp_2) == 0.9291776728412211

    @test RunStatistics.squares_cdf_approx(T_obs, L, epsp_2) == 0.9291776728412211
    @test RunStatistics.squares_pvalue_approx(T_obs, L, epsp_2) == 0.0708223271587789
end