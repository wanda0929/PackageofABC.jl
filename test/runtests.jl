using PackageofABC
using Yao
using TensorOperations
using LinearAlgebra

N = 4  # 系统大小
M = 3  # 生成算符的数量

# 定义 Pauli 矩阵
sigma_x = [0 1; 1 0]  # Pauli-X 矩阵
sigma_y = [0 -im; im 0]  # Pauli-Y 矩阵
sigma_z = [1 0; 0 -1]  # Pauli-Z 矩阵
I = [1 0; 0 1]  # 单位矩阵
# 定义 XXX 模型的 Hamiltonian
function xxx_hamiltonian(N::Int)
    H = zeros(Complex{Float64}, 2^N, 2^N)
    for i in 1:(N-1)
        H += kron(I^(i-1), sigma_x, sigma_x, I^(N-i-1)) +
             kron(I^(i-1), sigma_y, sigma_y, I^(N-i-1)) +
             kron(I^(i-1), sigma_z, sigma_z, I^(N-i-1))
    end
    return H
end

# 构建包含参数的生成算符
ops = [parameterized_operator(k[i]) for i in 1:M]

# 将生成算符表示为矩阵乘积算符
mpos = [mpo_from_operator(ops[i], N) for i in 1:M]

# 构建包含参数的量子线路
circuit = parameterized_circuit(k, N)

# 使用 QR 分解确保量子线路是幺正的
unitary_circuit = make_unitary(circuit)

# 构建 XXX 模型的 Hamiltonian
H = xxx_hamiltonian(N)

# 获取本征值和本征态
eigenvalues, eigenstates = get_eigenvalues_and_eigenstates(H)

# 打印本征值和本征态
println("Eigenvalues: ", eigenvalues)
println("Eigenstates: ", eigenstates)



