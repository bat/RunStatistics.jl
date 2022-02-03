# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).#

using RunStatistics
using Test
using Distributions


@testset "squares_additional" begin
    
    @testset "squares_paper" begin
        alpha = 0.001

        @test abs(squares_pvalue(15.5, 5) - alpha) <= 0.00002
        @test abs(squares_pvalue(23.8, 50) - alpha) <= 0.00002
        @test abs(squares_pvalue(25.6, 100) - alpha) <= 0.00008
    end

    @testset "squares_mathematica" begin
        n = 20
        T = [2, 5, 10, 20, 50]
        P = [0.8936721808595665, 0.42457437866154357, 0.06934906413527009, 0.0014159488909252227, 9.575188641974819e-9]
        eps = [1e-15, 1e-15, 1e-15, 1e-15, 1e-15]

        for i in 1:5
            @test abs(squares_pvalue(T[i], n) - P[i]) <= eps[i]
        end
    end

    @testset "squares_critical" begin
        N = 50 

        @test abs(squares_pvalue(19.645, 2 * N) - 0.01) <= 3e-5
        @test abs(squares_pvalue(15.34, 2 * N) - 0.05) <= 1e-4
    end
end

@testset "squares_approx_additional" begin
        
    @testset "squares_approx_split" begin

        Tobs = 15.5
        N = 12
        n = 2

        F = squares_cdf(Tobs, N)
        approx = F*F / (1.0 + RunStatistics.Delta(Tobs, N, N))

        @test abs(approx - squares_cdf_approx(Tobs, [N, n])) <= 1e-15
    end 

    function test_on_grid(K::Integer, N::Integer, n::Real, eps::Real)

        for i in 1:K
            Tobs = 20 + 2 * i
            @test abs(squares_cdf_approx(Tobs, [N, n]) - squares_cdf(Tobs, n * N)) <= eps
        end 
    end

    @testset "squares_approx_grid" begin

        test_on_grid(15, 40, 2, 2e-7)
        # test_on_grid(15, 40, 3, 4e-7) super expensive
    end 

    @testset "squares_approx_hH" begin
        
        Tobs = 15.5 
        x = 3.3
        N = 12 

        @test abs(RunStatistics.h(Tobs, N) - 0.000373964) <= 1e-8
        @test abs(RunStatistics.H(Tobs-x, Tobs, N) -  0.00245352) <= 1e-8
        @test abs(RunStatistics.Delta(Tobs, N, N) - 0.00175994) <= 1e-8
    end

    @testset "squares_approx_cdf" begin
        
        @test abs(cdf(Chisq(12), 15.5) - 0.784775) <= 1e-6
    end

    @testset "squares_approx_interpolate" begin
        Tobs = 32

        F71 = squares_cdf_approx(Tobs, [71, 5])
        F88 = squares_cdf_approx(Tobs, [88, 4])
        F89 = squares_cdf_approx(Tobs, [89, 4])
        F100 = squares_cdf_approx(Tobs, [Int64(100), 3.55])

        @test abs(F71 - F100) <= 3e-9
        @test abs(F71 - (F88 + 3*F89)/4) <= 2e-9
    end

    @testset "squares_approx_bound_error" begin

        Tobs = 8
        @test abs(squares_cdf(Tobs, 100) - (squares_cdf(Tobs, 50) ^ 2 - squares_cdf(Tobs, 100) * RunStatistics.Delta(Tobs, 100, 100))) <= 1e-3
    end
    
end
