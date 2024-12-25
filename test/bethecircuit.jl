using Test
using Yao
using LinearAlgebra 
using TensorOperations
using OMEinsum

# 导入 PackageofABC 模块
using PackageofABC

# Test the Pauli operators
@testset "pauli operators" begin
    # op_0 is magnon creation operator
    @test PackageofABC.op_0() == [1 0 0 0; 0 -1 1 0; 0 1 1 0; 0 0 0 -1]
end

#Test the contraction function
@testset "contraction" begin
    A = rand(2, 2, 2, 2)
    B = rand(2, 2, 2, 2)
    @test PackageofABC.contract(A, B, 1, 4) == tensorcontract(A, [0, -3, -2, -1], B, [1, 2, 3, 0])
end

# Test the contract2 function
@testset "contract2" begin
    C = rand(2, 2, 2)
    D = rand(2, 2, 2, 2)
    @test PackageofABC.contract2(C, D, 1, 4) == optein"aij,bkla->ijbkl"(C, D)
end

# Test the continued tensor contraction when M=1
@testset "tensor_from_operator" begin
    op0 = [1 0 0 0; 0 -1 1 0; 0 1 1 0; 0 0 0 -1]
    vars = rand(1)
    iden = [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    op1 = vars[1] * iden .+ im * op0
    @test PackageofABC.tensor_from_operator(1, vars) == reshape(op1, 2, 2, 2, 2)
end

# Test the continued tensor contraction when M=2
@testset "tensor_from_operator" begin
    op0 = [1 0 0 0; 0 -1 1 0; 0 1 1 0; 0 0 0 -1]
    vars = rand(2)
    iden = [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    op1 = vars[1] * iden .+ im * op0
    top1 = reshape(op1, 2, 2, 2, 2)
    op2 = vars[2] * iden .+ im * op0
    top2 = reshape(op2, 2, 2, 2, 2)
    final = PackageofABC.contract(top1, top2, 1, 4)
    @test PackageofABC.tensor_from_operator(2, vars) == final
end

# Test the continued tensor contraction when M=3
@testset "tensor_from_operator" begin
    op0 = [1 0 0 0; 0 -1 1 0; 0 1 1 0; 0 0 0 -1]
    vars = rand(3)
    iden = [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    op1 = vars[1] * iden .+ im * op0
    top1 = reshape(op1, 2, 2, 2, 2)
    op2 = vars[2] * iden .+ im * op0
    top2 = reshape(op2, 2, 2, 2, 2)
    final = PackageofABC.contract(top1, top2, 1, 4)
    op3 = vars[3] * iden .+ im * op0
    top3 = reshape(op3, 2, 2, 2, 2)
    final = PackageofABC.contract2(final, top3, 4, 4)
    @test PackageofABC.tensor_from_operator(3, vars) == final
end

# Test thetensor contraction when M=1
@testset "simplify_tensor" begin
    sim = rand(2, 2, 2, 2, 2, 2, 2, 2)
    @test PackageofABC.simplify_tensor(sim, 3) == sim[1, :, 1, 1, :, :, 1, :]
end

#Test the tensor contraction when M=2
@testset "simplify_tensor_2" begin
    sim = rand(2, 2, 2, 2, 2, 2, 2, 2)
    @test PackageofABC.simplify_tensor_2(sim, 3) == sim[:, :, 1, :, :, :, :, :]
end

# Test the tensor to matrix process
@testset "tensor_to_matrix2" begin
    tensor = rand(2, 2, 2, 2, 2, 2, 2)
    perm = permutedims(tensor, (4, 1, 2, 3, 5, 6, 7))
    @test PackageofABC.tensor_to_matrix2(tensor, (4, 1, 2, 3, 5, 6, 7)) == reshape(perm, 16, 8)
end

# Test the tensor to matrix process
@testset "tensor_to_matrix1" begin
    tensor = rand(2, 2, 2, 2, 2, 2)
    perm = permutedims(tensor, (4, 1, 2, 3, 5, 6))
    @test PackageofABC.tensor_to_matrix1(tensor, (4, 1, 2, 3, 5, 6)) == reshape(perm, 2, 32)
end

# test the unitary circuit
@testset "unitary_circuit" begin
    Mat1 = rand(2, 8)
    Mat2 = rand(16,8)
    vars = rand(3)
    R_1 = Mat1; #文中的G_0
    R_2 = Mat2; #第2个到第M个R_T
    M = length(vars)
    c = chain(M) #定义QR分解过程中的量子线路
    R_f = kron(eye(),R_1) * R_2 #将G0和R_T相乘
    Q, R = qr(R_f) #对获得的总矩阵进行QR分解
    P_i =Matrix(Q) #将获得的unitary矩阵提取出来
    R_1 = Matrix(R)  #将R矩阵提取出来，用于下一次迭代
    push!(c, put(M, (1,2)=>matblock(P_i))) #用unitary矩阵P_i构建量子线路
    R_f = kron(eye(),R_1) * R_2 #将G0和R_T相乘
    Q, R = qr(R_f) #对获得的总矩阵进行QR分解
    P_i = Matrix(Q)  #将获得的unitary矩阵提取出来
    R_1 = Matrix(R)  #将R矩阵提取出来，用于下一次迭代
    push!(c, put(M, (1,2,3)=>matblock(P_i))) #用unitary矩阵P_i构建量子线路
    p_f = transpose(Matrix(c))
    @test PackageofABC.unitary_circuit(Mat1, Mat2, vars) == p_f
end

# Test the Hamiltonian
@testset "hamiltonian" begin
    # Hamiltonian for 4 sites
    H = xxx_hamiltonian(4; periodic=false)
    # Check the Hamiltonian is Hermitian
    @test ishermitian(H)
    @test eigen(Matrix(xxx_hamiltonian(4; periodic=false))).values ≈ [-6.464101615137752, -3.828427124746195, -3.828427124746179,-3.8284271247461676, -1.0, -0.9999999999999991, -0.9999999999999987, 0.4641016151377566, 1.8284271247461898, 1.8284271247461907,1.8284271247461916, 3.0, 3.0, 3.0, 3.0, 3.0]
end

@testset "circuit" begin
    vars = rand(10)
    c = PackageofABC.circuit_from_operator(vars)
    @test length(c) == 3
    @test nqubits(c) == 4
    @info mat(c)
end

