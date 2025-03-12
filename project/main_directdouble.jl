using PackageofABC: circuitwithSingleblock, larger_circuit, circuit_from_operator
# 计算出参数值
# 定义多元方程组
#function equations!(F, u)
 #   F[1] = ((u[1] + im/2)/(u[1] - im/2))^N2 * ((u[1]-u[2]-im)/(u[1]-u[2]+im)) - 1
  #  F[2] = ((u[2] + im/2)/(u[2] - im/2))^N2 * ((u[2]-u[1]-im)/(u[2]-u[1]+im)) - 1
#end
N2 = 4
#function calculation_for_eigenvalue(N::Int)
    function equations!(F, u)
        F[1] = (u[1] + im/2)^N2 * (u[1] - u[1] - im) * (u[1] - u[2] - im) + (u[1] - im/2)^N2 * (u[1] - u[1] + im) * (u[1] - u[2] + im) 
        F[2] = (u[2] + im/2)^N2 * (u[2] - u[1] - im) * (u[2] - u[2] - im) + (u[2] - im/2)^N2 * (u[2] - u[1] + im) * (u[2] - u[2] + im) 
    end
    # 初始猜测值
    initial_guess = Complex{Float64}[1.0 + 0.0im, 1.0 + 0.0im]
    # 求解方程组
    result = nlsolve(equations!, initial_guess)    
    # 定义变量
    vars = result.zero
    M = length(vars)
    # 定义大块算符
    Matrices = circuit_from_operator(vars)
    # 定义量子线路
    circuit = larger_circuit(vars, NN, Matrices)
    # 初始化量子态为 |1000...⟩
    initial_amplitudes = zeros(Complex{Float64}, 2^(N2+M))
    initial_amplitudes[4] = 1.0  # |1000⟩ 对应的索引是 2 (从 1 开始计数)
    state = ArrayReg(initial_amplitudes)

    # 应用量子线路
    finalstate = apply!(state, circuit)

    # 测量量子态
    result = measure!(RemoveMeasured(), state, [1,2])

    normalize!(finalstate)
    realene = expect(xxx_hamiltonian(NN; periodic=true), finalstate)
    energy = 0
    for i in 1:M
        energy += 2/(4*vars[i]^2+1) 
    end
    if result == 0
        return energy - realene
    else
        return "Wrong"
    end
#end
calculation_for_eigenvalue(4)

0.5857864365030632 - 0.0im