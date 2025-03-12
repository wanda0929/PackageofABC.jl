using Test
using Yao
using LinearAlgebra 
using TensorOperations
using OMEinsum

# 导入 PackageofABC 模块
using PackageofABC
# Test the circuitwithSingleblock function
@testset "circuit" begin
   N = 4
   vars = rand(4)+im*rand(4)
    c = PackageofABC.circuitwithSingleblock(vars,N)
    @test length(c) == 4
    #@test nqubits(c) == 11
end

# Test the circuit_from_operator function
@testset "circuit" begin
    vars = rand(10) + im*rand(10)
    c = PackageofABC.circuit_from_operator(vars)
    @test length(c) == 4194304
    #@test nqubits(c) == 11
end

# test the larger_circuit function
@testset "larger_circuit" begin
    vars = rand(3) + im*rand(3)
    N = 4
    Matrices = rand(ComplexF64, 16,16)
    c = PackageofABC.larger_circuit(vars, N, Matrices)
    @test length(c) == 4
    #@test nqubits(c) == 7
end