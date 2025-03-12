using PackageofABC: circuitwithSingleblock, nlsolve, xxx_hamiltonian
using Yao
NN=2
#定义测量本征值的函数
function measure_eigenvalue(NN::Int)
    # 计算出参数值
    # 单magnon情形解方程
    function equation!(F, u)
        F[1] = (u[1] + im/2)^NN * (u[1] - im) - (u[1] - im/2)^NN * (u[1] + im)
    end
    # 初始猜测值
    initial_guess = Complex{Float64}[1.0 + 0.0im]
    # 求解方程组
    result = nlsolve(equation!, initial_guess)
    @show result
    # 定义变量vector
    vars = result.zero
    #vars = [10.0]
    # 定义量子线路
    circuit = circuitwithSingleblock(vars,NN)
    # 初始化量子态为 |1000⟩
    state = product_state(BitStr(vcat(zeros(Int, NN), [1])))
    # 应用量子线路
    apply!(state, circuit)
    # 测量ancella部分量子态
    result = measure!(RemoveMeasured(), state, 1)
    # 归一化量子态
    normalize!(state)
    # 计算哈密顿量下获得的量子态期望值
    caleig = expect(xxx_hamiltonian(NN; periodic=true), state)
    @show caleig
    # 数值计算XXXmodel本征能量
    simeig = 2/(4*vars[1]^2+1)-0.5
    result == 0 || @warn "Wrong"
    @show simeig caleig result
    return caleig - simeig
end

measure_eigenvalue(2)

    #0.535317324827776 - 1.017089829707777e-9im