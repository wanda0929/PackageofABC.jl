using PackageofABC: sigma_x, sigma_y, sigma_z, eye, zero, sigma_plus, sigma_minus, ident, contract, contract2, tensor_from_operator, xxx_hamiltonian, get_eigenvalues_and_eigenstates, tensor_to_matrix1,tensor_to_matrix2 , simplify_tensor, simplify_tensor_2, unitary_circuit 
using NLsolve
function energyofxxxmodel(N1::Int)
    
    # 定义多元方程组
    function equations!(F, u)
        #F[1] = ((u[1] + im/2)/(u[1] - im/2))^N1 * ((u[1]-u[2]-im)/(u[1]-u[2]+im)) * ((u[1]-u[3]-im)/(u[1]-u[3]+im)) * ((u[1]-u[4]-im)/(u[1]-u[4]+im)) - 1
        #F[2] = ((u[2] + im/2)/(u[2] - im/2))^N1 * ((u[2]-u[1]-im)/(u[2]-u[1]+im)) * ((u[2]-u[3]-im)/(u[2]-u[3]+im)) * ((u[2]-u[4]-im)/(u[2]-u[4]+im)) - 1
        #F[3] = ((u[3] + im/2)/(u[3] - im/2))^N1 * ((u[3]-u[1]-im)/(u[3]-u[1]+im)) * ((u[3]-u[2]-im)/(u[3]-u[2]+im)) * ((u[3]-u[4]-im)/(u[3]-u[4]+im)) - 1
        #F[4] = ((u[4] + im/2)/(u[4] - im/2))^N1 * ((u[4]-u[1]-im)/(u[4]-u[1]+im)) * ((u[4]-u[2]-im)/(u[4]-u[2]+im)) * ((u[4]-u[3]-im)/(u[4]-u[3]+im)) - 1
        F[1] = ((u[1] + im/2)^N1) * ((u[1]-u[1]-im) * (u[1]-u[2]-im) * (u[1]-u[3]-im) * (u[1]-u[4]-im)) + ((u[1] - im/2)^N1) * ((u[1]-u[1]+im) * (u[1]-u[2]+im) * (u[1]-u[3]+im) * (u[1]-u[4]+im))
        F[2] = ((u[2] + im/2)^N1) * ((u[2]-u[1]-im) * (u[2]-u[2]-im) * (u[2]-u[3]-im) * (u[2]-u[4]-im)) + ((u[2] - im/2)^N1) * ((u[2]-u[1]+im) * (u[2]-u[2]+im) * (u[2]-u[3]+im) * (u[2]-u[4]+im))
        F[3] = ((u[3] + im/2)^N1) * ((u[3]-u[1]-im) * (u[3]-u[2]-im) * (u[3]-u[3]-im) * (u[3]-u[4]-im)) + ((u[3] - im/2)^N1) * ((u[3]-u[1]+im) * (u[3]-u[2]+im) * (u[3]-u[3]+im) * (u[3]-u[4]+im))
        F[4] = ((u[4] + im/2)^N1) * ((u[4]-u[1]-im) * (u[4]-u[2]-im) * (u[4]-u[3]-im) * (u[4]-u[4]-im)) + ((u[4] - im/2)^N1) * ((u[4]-u[1]+im) * (u[4]-u[2]+im) * (u[4]-u[3]+im) * (u[4]-u[4]+im))
    end
    # 初始猜测值
    initial_guess = Complex{Float64}[1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im]
    # 求解方程组
    result = nlsolve(equations!, initial_guess)
    # 定义变量
    vars = result.zero
    M = length(vars)
    # 构建包含参数的生成算符
    ops = tensor_from_operator(M, vars)

    # 简化张量
    RT1 = simplify_tensor(ops, M)
    RT2 = simplify_tensor_2(ops, M)
    dim1 = size(RT1)
    dim2 = size(RT2)
    #张量转换成矩阵
    Mat1 = tensor_to_matrix1(RT1, Tuple(vcat(vcat(M,1:M-1),M+1)))
    Mat2 = tensor_to_matrix2(RT2,Tuple(vcat(vcat(2*M-1,1:2M-2),2*M, 2*M+1)) )

    #前N个张量连接的线路矩阵
    c_m = unitary_circuit(Mat1, Mat2, N1, vars)

    # 初始化量子态为 |1000⟩
    initial_amplitudes = zeros(Complex{Float64}, 2^(N1))
    initial_amplitudes[0b11110 + 1] = 1.0
    state = ArrayReg(initial_amplitudes)
    # 应用量子线路
    finalstate = apply!(state, c_m)
    # 归一化量子态
    normalize!(finalstate)
    # 计算哈密顿量下获得的量子态期望值
    caleig = expect(xxx_hamiltonian(N1; periodic=true), finalstate)
    # 数值计算XXXmodel本征能量
    simeig = -(2/(4*vars[1]^2+1) + 2/(4*vars[2]^2+1) + 2/(4*vars[3]^2+1) + 2/(4*vars[4]^2+1))
    return caleig - simeig
end
N1 = 5 # 系统大小
energyofxxxmodel(N1)



