# 定义 Pauli 矩阵
sigma_x() = [0 1; 1 0]  # Pauli-X 矩阵
sigma_y() = [0 -im; im 0]  # Pauli-Y 矩阵
sigma_z() = [1 0; 0 -1]  # Pauli-Z 矩阵
sigma_plus() = [0 1; 0 0]  # Pauli-升运算符
sigma_minus() = [0 0; 1 0]  # Pauli-降运算符
eye() = [1 0; 0 1]  # 单位矩阵
zero() = [0 0; 0 0]  # 零矩阵

# 构建包含参数的单块生成算符
op_0() = [sigma_z() sigma_minus(); sigma_plus() sigma_z()]
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
#iden = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]

function tensor_from_operator(M::Int, vars::Vector{Float64})
    op_1 = vars[1] * ident() .+ im * op_0()#单块含参量生成算符
    tensor_op_1 = reshape(op_1, 2, 2, 2, 2)#reshapeb-j-a-i
    if M == 1
        return tensor_op_1
    elseif M==2 
        op_2 = vars[2] * ident() .+ im * op_0()#第二块含参量生成算符
        tensor_op_2 = reshape(op_2, 2, 2, 2, 2)#reshapeb-j-a-i
        tensor_op_2 = contract(tensor_op_1, tensor_op_2, 1, 4)#收缩前两个张量
        return tensor_op_2
    else
        op_2 = vars[2] * ident() .+ im * op_0()
        tensor_op_2 = reshape(op_2, 2, 2, 2, 2)
        tensor_op_sum =  contract(tensor_op_1, tensor_op_2, 1, 4) #前i-1个的收缩结果
        for i in 3:M    
            op_i = vars[i] * ident() .+ im * op_0() #第i块含参量生成算符
            tensor_op_i = reshape(op_i , 2, 2, 2, 2) #把第i个reshapeb-j-a-i
            tensor_op_sum = contract2(tensor_op_sum, tensor_op_i, ((i-3)*2+4), 4)#将前i个的收缩结果代入第二个当中
            return tensor_op_sum
        end
        
    end
end
#test
vars = rand(3)
tensor_from_operator(3, vars)

# 简化第一个张量
#如果M=1，那么只有一个张量，不需要简化
function simplify_tensor(tensor, M::Int)
    input_index = 1# 固定第一个输入态为 0 态
    output_index = 1# 固定第一个输出态为 0 态
    vec = [output_index, :] # 用于重复的部分
    simplified_tensor = tensor[output_index, :, input_index, repeat(vec, M-2)..., :, output_index, :] #[j]-a-[i]-[j']-a'-[j'']-a''...b'''-[j''']-a'''
    return simplified_tensor #简化后的张量
end
# 简化第二个到第M个张量
function simplify_tensor_2(tensor, M::Int)
    input_index = 1# 固定第一个输入态为 0 态
    simplified_tensor = tensor[:, :, input_index, fill(:, 2*M-1)...] #j-a-[i]-j'-a'-j''-a''...b'''-j'''-a'''
    return simplified_tensor #简化后的张量
end
sim = rand(2, 2, 2, 2, 2, 2, 2, 2)
simplify_tensor_2(sim, 3)-sim[:, :, 1, :, :, :, :, :]


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

function xxx_hamiltonian(nbit::Int; periodic::Bool)
    map(1:(periodic ? nbit : nbit-1)) do i
        sum([kron(nbit, i=>operator,mod1(i+1, nbit)=>operator) for operator in [X, Y, Z]])
    end |> sum
end



