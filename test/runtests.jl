using LinearPowerFlow
using Test
using HiGHS
using JuMP


# # hack for local testing
# using Pkg
# Pkg.activate("..")
# using LinearPowerFlow
# Pkg.activate(".")


@testset "LinearPowerFlow.jl" begin
    

@testset "IEEE13 wye only" begin
    m = Model(HiGHS.Optimizer)
    p = Inputs(
        joinpath("data", "ieee13", "IEEE13Nodeckt.dss")
    )
    
end

end  # all tests
