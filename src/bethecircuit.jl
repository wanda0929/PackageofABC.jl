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

# 收缩相同维度的张量
function contract(A, B, kA, kB)
    # 获取张量 A 的维度索引，范围从 -ndims(A) 到 -1
    IA = collect(-ndims(A):-1)
    
    # 获取张量 B 的维度索引，范围从 1 到 ndims(B)
    IB = collect(1:ndims(B))
    
    # 将张量 A 的第 kA 个维度索引设置为 0
    IA[kA] = 0
    
    # 将张量 B 的第 kB 个维度索引设置为 0
    IB[kB] = 0
    
    # 使用 tensorcontract 函数对张量 A 和 B 进行收缩
    return tensorcontract(A, IA, B, IB)
end

# 收缩不同维度的张量
function contract2(C,D,kC,kD)
    IC = ntuple(d -> d==kC ? 0 : -d, ndims(C))
    ID = ntuple(d -> d==kD ? 0 : +d, ndims(D))
    return TensorOperations.tensorcontract(C,IC,D,ID)
end

#将构建的生成算符连接成张量
function tensor_from_operator(M::Int, vars::Vector{Complex{Float64}})
    op_1 = vars[1] * ident() .+ im * op_0()#单块含参量生成算符
    tensor_op_1 = reshape(op_1, 2, 2, 2, 2)#reshapeb-j-a-i
    if M == 1
        return tensor_op_1
    elseif M==2 
        op_2 = vars[2] * ident() .+ im * op_0()#第二块含参量生成算符
        tensor_op_2 = reshape(op_2, 2, 2, 2, 2)#reshapeb-j-a-i
        tensor_op_2 = contract(tensor_op_1, tensor_op_2, 2, 4)#收缩前两个张量
        return tensor_op_2
    else
        op_2 = vars[2] * ident() .+ im * op_0()
        tensor_op_2 = reshape(op_2, 2, 2, 2, 2)
        tensor_op_sum =  contract(tensor_op_1, tensor_op_2, 2, 4) #前i-1个的收缩结果
        for i in 3:M    
            op_i = vars[i] * ident() .+ im * op_0() #第i块含参量生成算符
            tensor_op_i = reshape(op_i , 2, 2, 2, 2) #把第i个reshapeb-j-a-i
            tensor_op_sum = contract2(tensor_op_sum, tensor_op_i, (2*i-1), 4)#将前i个的收缩结果代入第二个当中
        end  
        return tensor_op_sum
    end
end

Rmatrix(var) = var * ident() .+ im * op_0()
function circuit_from_operator(vars::Vector{Complex{Float64}})
    M = length(vars)
    c = chain(M+1)
    for i in 1:M
        # `Rmatrix` returns a 4x4 matrix
        # `matblock` is used to convert it to a gate on 2 qubits
        # `put` is used to elivate it to a gate on M qubits
        push!(c, put(M+1, (i, i+1)=>matblock(Rmatrix(vars[i]))))
    end
    return c
end
#test
#vars = rand(3)
#tensor_from_operator(3,vars)




# 简化第一个张量
#如果M=1，那么只有一个张量，不需要简化
function simplify_tensor(tensor, M::Int)
    input_index = 1# 固定第一个输入态为 0 态
    output_index = 1# 固定第一个输出态为 0 态
    vec = [output_index, :] # 用于重复的部分
    simplified_tensor = tensor[output_index, :, input_index, repeat(vec, M-2)..., output_index, :, :] #[b]-a-[i]-[b']-a'-[b'']-a''...b'''-[j''']-a'''
    return simplified_tensor #简化后的张量
end

# 简化第二个到第M个张量
function simplify_tensor_2(tensor, M::Int)
    input_index = 1# 固定第一个输入态为 0 态
    simplified_tensor = tensor[:, :, input_index, fill(:, 2*M-1)...] #b-a-[i]-b'-a'-b''-a''...b'''-j'''-a'''
    return simplified_tensor #简化后的张量
end

# 将高维张量R1变回矩阵
function tensor_to_matrix1(tensor, new_dims::NTuple{N, Int}) where {N}
    # 重新排列张量的维度，使得要合并的维度在一起
    tensor_permuted = permutedims(tensor, new_dims)
    # 获取张量的尺寸
    dims = size(tensor_permuted)
    # 计算合并后的维度
    merged_dim = prod(dims[1:N-1])
    single_dim_size = dims[N]
    
    # 将张量重新排列成矩阵
    return reshape(tensor_permuted, single_dim_size, merged_dim)
end

# 将高维张量R2-RM变回矩阵
function tensor_to_matrix2(tensor, new_dims::NTuple{N, Int}) where {N}
    # 重新排列张量的维度，使得要合并的维度在一起
    tensor_permuted = permutedims(tensor, new_dims)
    # 获取张量的尺寸
    dims = size(tensor_permuted)
    # 计算合并后的维度
    merged_dim = prod(dims[1:(N+1)/2])
    single_dim =prod(dims[(N+3)/2:N])
    
    # 将张量重新排列成矩阵
    return reshape(tensor_permuted, merged_dim, single_dim)
end
#tensor = rand(2, 2, 2, 2, 2, 2, 2)
#tensor_to_matrix2(tensor, (4, 1, 2, 3, 5, 6, 7))



# 生成unitary矩阵的量子线路
function unitary_circuit(Matrix1, Matrix2, N::Int, vars::Vector{Complex{Float64}})
    R_1 = Matrix1; #文中的G_0
    R_2 = Matrix2; #第2个到第N个R_T
    M = length(vars)
    Matrices = Dict()
    c = chain(N) #定义QR分解过程中的量子线路
    for i in 1:N-1
        R_f = kron(eye(),R_1) * R_2 #将G0和R_T相乘
        Q, R = qr(R_f) #对获得的总矩阵进行QR分解
        P_i = Matrix(Q) #将获得的unitary矩阵提取出来
        R_1 = Matrix(R) #将R矩阵提取出来，用于下一次迭代
        Matrices["matblock_$i"] = P_i
        #push!(c, put(M, vcat(1:i+1)=>matblock(P_i))) #用unitary矩阵P_i构建量子线路
    end

    for i in 1:N-M
        push!(c, put(N, vcat(N-M+2-i:N+1-i)=>matblock(Matrices["matblock_$(N-i)"]))) #用unitary矩阵P_i构建量子线路
    end

    for i in N-M+1:N-1
        push!(c, put(N, vcat(1:N-i+1)=>matblock(Matrices["matblock_$(N-i)"])))
         #用unitary矩阵P_i构建量子线路
    end
    return c #返回前M-1个unitary矩阵组成的量子线路
end


function lager_unitary_circuit(Matrix1, Matrix2, N::Int, vars::Vector{Complex{Float64}})
    R_1 = Matrix1; #文中的G_0
    R_2 = Matrix2; #第2个到第N个R_T
    M = length(vars)
    Matrices = Dict()
    c = chain(N) #定义QR分解过程中的量子线路
    for i in 1:N-1
        R_f = kron(eye(),R_1) * R_2 #将G0和R_T相乘
        Q, R = qr(R_f) #对获得的总矩阵进行QR分解
        P_i = Matrix(Q) #将获得的unitary矩阵提取出来
        R_1 = Matrix(R) #将R矩阵提取出来，用于下一次迭代
        Matrices["matblock_$i"] = P_i
        #push!(c, put(M, vcat(1:i+1)=>matblock(P_i))) #用unitary矩阵P_i构建量子线路
    end
    for i in 1:N-M+1
        push!(c, put(M, vcat((N-M)+2-i:N-i+1)=>matblock(Matrices["matblock_$(N-i)"]))) #用unitary矩阵P_i构建量子线路
    end

    for i in N-M+1:N
        push!(c, put(M, vcat(1:N+1-i)=>matblock(Matrices["matblock_$(N-i)"]))) #用unitary矩阵P_i构建量子线路
    end
    return c #返回前M-1个unitary矩阵组成的量子线路
end

#


# 获取 XXX 模型的本征值和本征态
function get_eigenvalues_and_eigenstates(H)
    eigenvalues, eigenstates = eigen(H)
    return eigenvalues, eigenstates
end

# xxx模型哈密顿量
function xxx_hamiltonian(nbit::Int; periodic::Bool)
    map(1:(periodic ? nbit : nbit-1)) do i
        sum([kron(nbit, i=>operator,mod1(i+1, nbit)=>operator) for operator in [X, Y, Z]])
    end |> sum
end

#test