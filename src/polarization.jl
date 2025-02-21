using LinearAlgebra
using WignerD

"""
ElectricVector(pol, θ, ϕ)

DESCRIPTION
Computes the light vector E from Cartesian to Cyclic coordinates, to be compatable with Wigner theorem.

INPUT
E vector in Cartesian coordinates
pol: polarization of the light - +1, 0 or -1
θ: rotation angle about y-axis
ϕ: rotation angle about z'-axis

OUTPUT
E vector in Cyclic coordinates with 3 q components: -1, 0 and +1 respectively
"""

struct ElectricVector

    pol::Vector{ComplexF64}
    α::Float64
    β::Float64
    γ::Float64
    cyclic::Array{Complex{Float64},1}

    function ElectricVector(pol, α, β, γ)
    

        cyclic = wignerD(1,α,β,γ)' * pol

        new(pol, α, β, γ, cyclic)
    end
end