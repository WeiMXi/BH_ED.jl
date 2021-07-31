using BH_ED
using Test

@testset "BH_ED.jl" begin
    r1 = BH_solver(N   = 6
                , Lp   = 3
                , L    = 3*6
                , μ    = 1.0
                , J    = 1.0
                , V    = 1.0
                , ε    = 1.0
    )
    @test abs(r1[1][1] + 8.23423318868652) <= 0.1
end
