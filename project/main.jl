using PackageofABC: sigma_x, sigma_y, sigma_z, eye, parameterized_operator, parameterized_circuit, make_unitary, get_eigenvalues_and_eigenstates
N = 4  # 系统大小
M = 3  # 生成算符的数量

using Yao

using Test, LinearAlgebra
# 构建包含参数的生成算符
ops = [parameterized_operator(k[i]) for i in 1:M]

# 将生成算符表示为矩阵乘积算符
mpos = [mpo_from_operator(ops[i], N) for i in 1:M]

#...

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
