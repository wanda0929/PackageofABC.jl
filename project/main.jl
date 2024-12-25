using PackageofABC: sigma_x, sigma_y, sigma_z, eye, zero, sigma_plus, sigma_minus, ident, contract, contract2, tensor_from_operator, xxx_hamiltonian, get_eigenvalues_and_eigenstates, tensor_to_matrix1,tensor_to_matrix2 , simplify_tensor, simplify_tensor_2
N = 4  # 系统大小
M = 3  # 生成算符的数量
# 构建包含参数的生成算符
ops = tensor_from_operator(M, rand(M))

# 简化张量
RT1 = simplify_tensor(ops, M)
RT2 = simplify_tensor_2(ops, M)
dim1 = size(RT1)
dim2 = size(RT2)
#张量转换成矩阵
Mat1 = tensor_to_matrix1(RT1, Tuple(vcat(vcat(M,1:M-1),M+1)))
Mat2 = tensor_to_matrix2(RT2,Tuple(vcat(vcat(2*M-1,1:2M-2),2*M, 2*M+1)) )

#QR分解






g0Rt = kron(ident(),Mat1)
Mat2
Q, R = qr(g0Rt)
for i in 1:M
    Q, R = qr(Mat1[:,:,i])
    Mat1[:,:,i] = Q
end

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
