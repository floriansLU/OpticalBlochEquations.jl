using Distributed
@everywhere begin
    
using OpticalBlochEquations
using Test
using Plots
gr()
using DelimitedFiles


#include("../src/param.jl")
#include("../src/polarization.jl")
#include("../src/AtomicTransitions.jl")


@testset "OpticalBlochEquations.jl" begin
    # Write your tests here.

    params = param(cezijsD1)
    laser = lazers()
    #n2Fm_ats_g = fill_n2Fm_ats_g(par)
    #n2Fm_ats_e = fill_n2Fm_ats_e(par)

    #gDict = fill_gDict(n2Fm_ats_g)
    #eDict = fill_eDict(n2Fm_ats_e)



#ierosmes, novērošanas un zondēšanas ģeometrijas definēšana (pol, θ, ϕ)
    e_vec_i = ElectricVector(1, π / 2, 0).cyclic
    e_vec_n = ElectricVector(1, π / 2, π / 2).cyclic
    e_vect_z = ElectricVector(0, π / 2, π / 4).cyclic
    #dip_star = lin_fillDipoleMatrix_star(par, e_vec_i, n2Fm_ats_g, n2Fm_ats_e) #izveido dipolu matricu
    #dip = lin_fillDipoleMatrix(par, e_vec_i, n2Fm_ats_g, n2Fm_ats_e)
    #GAMMA = fill_GAMMA(par, laz, n2Fm_ats_g, n2Fm_ats_e) #izveido spontānās relaksācijas matricu

#induced relaxation
#    λ = zeros(Complex{Float64}, par.dim_g * par.dim_g + par.dim_e * par.dim_e)
#    for index = axes(gDict, 1)
#        gi = gDict[index, 1]
#        gj = gDict[index, 2]
#        if (gi == gj)
#            λ[index] = laz.γ / par.dim_g
#        end
#    end

"""Dopler shift"""
    #step_lit = 2
    #sigma = sqrt(laz.kB * laz.T / laz.masa) * laz.ω_svitr / laz.c
    #dscan = 2 * sigma * step_lit
    #dstep_nr = 150 * step_lit #number of steps
    #dstep_length = dscan / dstep_nr

    #Brange = [-40:2:-1; -0.9:0.1:0] |> mirror_around_end
    #Brange=[-40,-30,-20,-10,-5,-1,-0.1,0,0.1,5,10,20, 30, 40]
    #Brange=[-40,-20,-5,-1,0,1,5,20,40]
    #Brange=[-2000,-1500,-1000,-500,-100,-40,-1,-0.1,0,0.1,1,40,100,500,1000,1500,2000]
    #Brange=[-1,0,1]
    #res = @timed pmap(signals_for_pmap, Brange)   ### FHG 2024-02-16
    B₀=0
    evecs=(e_vec_i,e_vec_n,e_vect_z)
    Doppler_steps=150
    #signals(B₀, par, laz, gDict, eDict, n2Fm_ats_g, n2Fm_ats_e)
    #signals(B₀, params, laser, evecs, Doppler_steps)
    signals_for_pmap(B₀) = signals(B₀, params, laser, evecs, Doppler_steps)
    #Brange=[-5.0,0.0,5.0]
    Brange=[-40,-35,-30,-25,-20,-15,-10,-8,-6,-5,-4,-3,-2,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,
    0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,8,10,15,20,25,40]
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
    savefig(p, "A-ar-doplera-efektu.png")

    p = plot(Brange, CsI)
    plot!(p, title="I ar Doplera efektu", xlabel="B")
    savefig(p, "I-ar-doplera-efektu.png")

    open("I_Doplera.txt", "w") do io       ### FHG added 2024-02-02
        writedlm(io, [Brange CsI] )               ### FHG added 2024-02-02
    end                                    ### FHG added 2024-02-02

    open("A_Doplera.txt", "w") do io       ### FHG added 2024-02-02
        writedlm(io, [Brange CsA] )                ### FHG added 2024-02-02
    end                                    ### FHG added 2024-02-02

end # end @test

end # end @everywhere
