module PackageofABC

using LinearAlgebra
using TensorOperations

# 定义 Pauli 矩阵
sigma_x = [0 1; 1 0]  # Pauli-X 矩阵
sigma_y = [0 -im; im 0]  # Pauli-Y 矩阵
sigma_z = [1 0; 0 -1]  # Pauli-Z 矩阵
sigma_plus = [0 1; 0 0]  # Pauli-升运算符
sigma_minus = [0 0; 1 0]  # Pauli-降运算符
I = [1 0; 0 1]  # 单位矩阵


# 构建包含参数的生成算符
function parameterized_operator(M::Int, vars::Vector{Float64})
    op_0 = [sigma_z sigma_minus; sigma_plus sigma_z]
    op = I
    for i in 1:M
        op_1 = vars[i] + im * op_0
        op *= kron(I^(M+1-i), op_1, I^(i-1))
    end
    return op
end

# 将生成算符表示为张量
function tensor_from_operator(op::Matrix{Complex{Float64}}, M::Int)
    d = 2^(M+1)
    tensor_op = reshape(op, ntuple(i->2, 2M+2)...)
    return tensor_op
end

# 简化张量
# 固定输入态和输出态
function simplify_tensor(tensor, M::Int)
    # 固定第一个输入态为 0 态
    input_index = 1
    output_index = 1
    simplified_tensor = tensor[input_index, output_index...]
    # 生成从 i2 到 i(M+1) 的索引
    input_indices = Symbol.("i" .* string.(1:M+1))
    input_diff = Symbol.("i" .* string.(1))
    output_indices = Symbol.("j" .* string.(1:M+1))
    new_input_indices = Symbol.("i" .* string.(2:M+1))

    # 简化张量
    simplified_dims = ntuple(_ -> 2, M + 1)
    simplified_tensor = Array{Complex{Float64}, M+1}(undef, simplified_dims...)
    
    @tensor begin
        simplified_tensor[new_input_indices..., output_indices...] := tensor[input_indices..., output_indices...] * input_state[input_diff]
        for k in 1:M
            new_output_indices = output_indices[1+k:M+1]
            simplified_tensor[input_indices..., new_output_indices] = simplified_tensor[input_indices..., output_indices[k], new_output_indices]*output_state[output_indices[k]]
        end
    end
    
    return simplified_tensor
end

# 将高维张量变回矩阵
function tensor_to_matrix(tensor, input_dims::NTuple{N, Int}, output_dims::NTuple{M, Int}) where {N, M}
    input_size = prod(input_dims)
    output_size = prod(output_dims)
    return reshape(tensor, input_size, output_size)
end

# 构建包含参数的量子线路
function parameterized_circuit(k::Vector{Float64}, N::Int)
    circuit = chain(N)
    
    
    return circuit
end

# 使用 QR 分解确保量子线路是幺正的
function make_unitary(circuit)
    U, R = qr(mat(circuit))
    unitary_circuit = matblock(U)
    return unitary_circuit
end

# 获取 XXX 模型的本征值和本征态
function get_eigenvalues_and_eigenstates(H)
    eigenvalues, eigenstates = eigen(H)
    return eigenvalues, eigenstates
end










end
