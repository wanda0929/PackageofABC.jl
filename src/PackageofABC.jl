module PackageofABC

using OMEinsum
using LinearAlgebra
using Yao
using TensorOperations
using NLsolve
export op_0,eye, ident, contract, contract2, tensor_from_operator, simplify_tensor, simplify_tensor_2, tensor_to_matrix, xxx_hamiltonian,get_eigenvalues_and_eigenstates, unitary_circuit

include("bethecircuit.jl")
include("Directbethe.jl")

end