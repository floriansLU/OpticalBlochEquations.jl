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
    pol::Int
    θ::Float64
    ϕ::Float64
    #α::Float64
    #β::Float64
    #γ::Float64
    cyclic::Array{Complex{Float64},1}

    function ElectricVector(pol, θ, ϕ)
    #function ElectricVector(pol, α, β, γ)
        #Transformāciju m-ca starp Dekarta un sfērisko bazi
        U = [-(1 / sqrt(2)) im/sqrt(2) 0; 0 0 1; 1/sqrt(2) im/sqrt(2) 0]

        # if (pol == 1) # Pa kreisi polarizēts 1/sqrt(2)(x + iy)
        #     initial = [1 / sqrt(2), (1 / sqrt(2))im, 0]
        # elseif (pol == 0) # Lineaars z
        #     initial = [0, 0, 1]
        # elseif (pol == -1) # Pa labi polarizēts 1/sqrt(2)(x - iy)
        #     initial = [1 / sqrt(2), -(1 / sqrt(2))im, 0]
        # end

        initials = Dict(1  => [-1 / sqrt(2), (-1 / sqrt(2))im, 0], # FHG 2024-07-22: added minus signs
                        0  => [0, 0, 1],
                        -1 => [1 / sqrt(2), -(1 / sqrt(2))im, 0])

        ## contravariant spherical vectors
        #initials = Dict(1  => [-1 / sqrt(2), (1 / sqrt(2))im, 0], # FHG 2024-07-22: added minus signs
        #                0  => [0, 0, 1],
        #                -1 => [1 / sqrt(2), (1 / sqrt(2))im, 0])

        #pagrieziena m-ca ap y asi 
        R1 = [cos(θ) 0 sin(θ); 0 1 0; -sin(θ) 0 cos(θ)]

        #pagrieziena m-ca ap z asi
        R2 = [cos(ϕ) -sin(ϕ) 0; sin(ϕ) cos(ϕ) 0; 0 0 1]

        #Kopējā transformācijas matrica
        R = R2 * R1

        #Dekarta vektora transformācija
        rotated_cart = R * initials[pol]

        #Pareja no Dekarta uz cikliskajam koordinātām
        cyclic = U * rotated_cart
        #cyclic = initials[pol] * wignerD(1,α,β,γ)

        new(pol, θ, ϕ, cyclic)
        #new(pol, α, β, γ, cyclic)
    end
end