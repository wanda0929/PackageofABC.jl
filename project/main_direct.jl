using PackageofABC: circuitwithSingleblock, larger_circuit, circuit_from_operator, xxx_hamiltonian
using NLsolve
using Yao
# 计算出参数值
N2=5
function calculation_for_eig(N2::Int)
    # 定义多元方程组
    
    function equations!(F, u)
        F[1] = ((u[1] + im/2)/(u[1] - im/2))^N2 * ((u[1]-u[2]-im)/(u[1]-u[2]+im)) * ((u[1]-u[3]-im)/(u[1]-u[3]+im)) * ((u[1]-u[4]-im)/(u[1]-u[4]+im)) - 1
        F[2] = ((u[2] + im/2)/(u[2] - im/2))^N2 * ((u[2]-u[1]-im)/(u[2]-u[1]+im)) * ((u[2]-u[3]-im)/(u[2]-u[3]+im)) * ((u[2]-u[4]-im)/(u[2]-u[4]+im)) - 1
        F[3] = ((u[3] + im/2)/(u[3] - im/2))^N2 * ((u[3]-u[1]-im)/(u[3]-u[1]+im)) * ((u[3]-u[2]-im)/(u[3]-u[2]+im)) * ((u[3]-u[4]-im)/(u[3]-u[4]+im)) - 1
        F[4] = ((u[4] + im/2)/(u[4] - im/2))^N2 * ((u[4]-u[1]-im)/(u[4]-u[1]+im)) * ((u[4]-u[2]-im)/(u[4]-u[2]+im)) * ((u[4]-u[3]-im)/(u[4]-u[3]+im)) - 1
    end
    # 初始猜测值
    initial_guess = Complex{Float64}[1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im]
    # 求解方程组
    result = nlsolve(equations!, initial_guess)
    # 定义变量
    vars = result.zero
    M = length(vars)
    # 定义大块算符
    Matrices = circuit_from_operator(vars, N2)
    # 定义量子线路
    circuit = larger_circuit(vars, N2, Matrices)
    # 初始化量子态为 |1000...⟩
    initial_amplitudes = zeros(Complex{Float64}, 2^(N2+M))
    initial_amplitudes[0b111100000] = 1.0  # |1000⟩ 对应的索引是 2 (从 1 开始计数)
    state = ArrayReg(initial_amplitudes)
    # 应用量子线路
    finalstate = apply!(state, circuit)
    # 测量量子态
    result = measure!(RemoveMeasured(), state, [1,2,3,4])
    normalize!(finalstate)
    realene = expect(xxx_hamiltonian(N2; periodic=true), finalstate)
    energy = 0
    for i in 1:M
        energy += 2/(4*vars[i]^2+1) 
    end
    energy
    if result == 0000
        return realene
    else
        return "Wrong"
    end
end
#calculation_for_eig(5)
while true
    output = calculation_for_eig(N2)
    if typeof(output) != String
        println("Result: ", output)
        break
    end
end