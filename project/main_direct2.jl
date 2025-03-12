using PackageofABC: circuitwithSingleblock, larger_circuit, circuit_from_operator
N2 = 4
function calculation_for_eigenvalue(N::Int)
    vars = [1/(2*sqrt(3)), -1/(2*sqrt(3))]
    M = length(vars)
    # 定义大块算符
    Matrices = circuit_from_operator(vars)
    # 定义量子线路
    circuit = larger_circuit(vars, N2, Matrices)
    # 初始化量子态
    #state = product_state(BitStr(vcat(zeros(Int, N2), ones(Int, M))))
    initial_amplitudes = zeros(Complex{Float64}, 2^(N2+M))
    initial_amplitudes[0b000011+1] = 1.0  # |1000⟩ 对应的索引是 2 (从 1 开始计数)
    state = ArrayReg(initial_amplitudes)

    # 应用量子线路
    finalstate = apply!(state, circuit)

    # 测量量子态
    result = measure!(RemoveMeasured(), state, [1,2])
    @show result

    normalize!(finalstate)
    realene = expect(xxx_hamiltonian(N2; periodic=true), finalstate)
    energy = 0
    for i in 1:M
        energy += 2/(4*vars[i]^2+1) 
    end
    energy = energy - 1
    if result == 0
        return energy - realene
    else
        return "Wrong"
    end
end
calculation_for_eigenvalue(4)




