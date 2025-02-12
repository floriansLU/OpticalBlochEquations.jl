"""Test Probability Surfaces"""

using OpticalBlochEquations
using Test
using Plots
gr()

#include("../src/param.jl")
#include("../src/polarization.jl")
#include("../src/AtomicTransitions.jl")
#include("../src/ProbabilitySurfaces.jl")

@testset "OpticalBlochEquations.jl" begin
    # Write your tests here.

    params = param(cesiumD1);
    laser_params = laser();




#ierosmes, novērošanas un zondēšanas ģeometrijas definēšana (pol, θ, ϕ)
#    e_vec_ex = ElectricVector(1, π / 2, 0).cyclic
#    e_vec_obs = ElectricVector(1, π / 2, π / 2).cyclic
#    e_vec_probe = ElectricVector(0, π / 2, π / 4).cyclic
    e_vec_ex = ElectricVector(1, 0, π / 2, 0).cyclic
    e_vec_obs = ElectricVector(1, 0, π / 2, π / 2).cyclic
    e_vec_probe = ElectricVector(0, 0, π / 2, π / 4).cyclic





    B₀=1000
    evecs=(e_vec_ex,e_vec_obs,e_vec_probe)
    Doppler_steps=150
    #signals(B₀, par, laz, gDict, eDict, n2Fm_ats_g, n2Fm_ats_e)
    #signals(B₀, params, laser, evecs, Doppler_steps)
    signals_for_pmap(B₀) = signals(B₀, params, laser_params, evecs, Doppler_steps)
    #Brange=[-5.0,0.0,5.0]
    res=signals_for_pmap(B₀)
    ρgg=res[3]
    ρee=res[4]
    ρggJ4=ρgg[1:9,1:9]
    J=4
    surface=plotProbSurf(J,ρggJ4)
    show(surface)
end