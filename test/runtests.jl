# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).

import Test
import RunStatistics
import Documenter

Test.@testset "Package RunStatistics" begin
    include("test_Partitions.jl")
    include("test_Tobs.jl")
    include("test_Squares.jl")
    # doctests
    Documenter.DocMeta.setdocmeta!(
        RunStatistics,
        :DocTestSetup,
        :(using RunStatistics);
        recursive = true
    )
    Documenter.doctest(RunStatistics)
end # testset


