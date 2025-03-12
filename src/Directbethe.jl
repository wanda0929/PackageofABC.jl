# 定义 Pauli 矩阵
sigma_x() = [0 1; 1 0]  # Pauli-X 矩阵
sigma_y() = [0 -im; im 0]  # Pauli-Y 矩阵
sigma_z() = [1 0; 0 -1]  # Pauli-Z 矩阵
sigma_plus() = [0 1; 0 0]  # Pauli-升运算符
sigma_minus() = [0 0; 1 0]  # Pauli-降运算符
eye() = [1 0; 0 1]  # 单位矩阵
zero() = [0 0; 0 0]  # 零矩阵

# 构建包含参数的单块生成算符
op_0() = [sigma_z() sigma_minus(); sigma_plus() -sigma_z()]
ident() = [eye() zero(); zero() eye()]
#op_0 = [1 0 0 0; 0 -1 1 0; 0 1 1 0; 0 0 0 -1]
Rmatrix(var) = var * ident() .+ im * op_0()

# The circuit for onelayer
function circuitwithSingleblock(vars::Vector{Float64}, N)
    c = chain(N+1)
    for i in 1:N
        # `Rmatrix` returns a 4x4 matrix
        # `matblock` is used to convert it to a gate on 2 qubits
        # `put` is used to elivate it to a gate on M qubits
        push!(c, put(N+1, (i, i+1)=>matblock(Rmatrix(vars[1]))))
    end
    return c
end




function circuit_from_operator(vars::Vector{Float64})
    M = length(vars)
    c = chain(M+1)
    for i in 1:M
        # `Rmatrix` returns a 4x4 matrix
        # `matblock` is used to convert it to a gate on 2 qubits
        # `put` is used to elivate it to a gate on M qubits
        push!(c, put(M+1, (i, i+1)=>matblock(Rmatrix(vars[i]))))
    end
    Matrices = Matrix(c)
    return Matrices
    
end

# 进一步构建circuit
function larger_circuit(vars::Vector{Float64}, N::Int, Matrices)
    M = length(vars)
    c = chain(N+M)
    for i in 1:N
        push!(c, put(M+N, Tuple(vcat(N+1-i: M+N+1-i))=>matblock(Matrices)))
    end
    return c
end


