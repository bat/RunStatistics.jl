# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).

using RunStatistics
using Test

@testset "t_obs" begin

    @test RunStatistics.t_obs([-1, 1, 3, -2], 0, 1) == (10.0, [2,3])
    @test RunStatistics.t_obs([-1, 1, 3, -2], zeros(4), ones(4)) == (10.0, [2,3])

    @test RunStatistics.t_obs([-1, 1, 3, -2, 1, 3], 0, 1) == (10.0, [[2,3],[5, 6]])
    @test RunStatistics.t_obs([-1, 1, 3, -2, 1, 3], zeros(6), ones(6)) == (10.0, [[2,3],[5, 6]])
end