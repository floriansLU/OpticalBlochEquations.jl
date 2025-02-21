using Distributed
@everywhere begin
    
using OpticalBlochEquations
using Test
using Plots
gr()
using DelimitedFiles





@testset "OpticalBlochEquations.jl" begin
    # Write your tests here.

    params = param(cesiumD1)
    laser_params = laser()




# Definition of excitation, fluorescence and probe polarization using spherical coordinates and Euler angles. 

    e_vec_ex = ElectricVector([1,0,0], 0, π / 2, 0).cyclic
    e_vec_obs = ElectricVector([1,0,0], 0, π / 2, π / 2).cyclic
    e_vec_probe = ElectricVector([0,0,0], 0, π / 2, π / 4).cyclic
 

"""Doppler shift"""
    
    B₀=0
    evecs=(e_vec_ex,e_vec_obs,e_vec_probe)
    Doppler_steps=150
    signals_for_pmap(B₀) = signals(B₀, params, laser_params, evecs, Doppler_steps)
    Brange=[-5.0,0.0,5.0]
    #Brange=[-40,-35,-30,-25,-20,-15,-10,-8,-6,-5,-4,-3,-2,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,
    #0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,8,10,15,20,25,40]
    res = @timed pmap(signals_for_pmap, Brange)
    CsI=[]
    CsA=[]

    for i in 1:length(Brange)
        push!(CsA,res[1][i][2] |> real)
        push!(CsI,res[1][i][1] |> real)
    end

    default(legend=false)
    p = plot(Brange, CsA)
    plot!(p, title="A ar Doplera efektu", xlabel="B")
    #savefig(p, "A-ar-doplera-efektu.png")

    p = plot(Brange, CsI)
    plot!(p, title="I ar Doplera efektu", xlabel="B")
    #savefig(p, "I-ar-doplera-efektu.png")

  

end # end @test

end # end @everywhere
