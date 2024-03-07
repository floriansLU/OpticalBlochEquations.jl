module OpticalBlochEquations

export signals, signals_for_pmap, fill_n2Fm_ats_e, fill_n2Fm_ats_g, param, lazers

# Write your package code here.

using SparseArrays
using QuantumOptics
using WignerSymbols
using QuantumOptics.steadystate
using Plots
gr()
using DelimitedFiles
using BenchmarkTools
using Distributed

include("../src/param.jl")
include("../src/polarization.jl")

""" Lande Factor """
#function LandeFactorJ(J,L,S)
#    return 1+(J*(J+1)+S*(S+1)-L*(L+1))/(2*J*(J+1))
#end
#

#function LandeFactorF(F,I,J, gJ) 
#    return gJ*(F*(F+1)-I*(I+1)+J*(J+1))/(2*F*(F+1))
#end

function fill_n2Fm_ats_g(par::param)
    Fg_array = Int[]
    mg_array = Int[]
    for Fg = Int(par.grFmax):-1:Int(par.grFmin)
        Fg_array = [Fg_array; fill(Fg, 2*Fg + 1)]
        mg_array = [mg_array; Fg:-1:-Fg]
    end
    return hcat(Fg_array, mg_array)
end

function fill_n2Fm_ats_e(par::param)
    Fe_array = Int[]
    me_array = Int[]
    for Fe = Int(par.exFmax):-1:Int(par.exFmin)
        Fe_array = [Fe_array; fill(Fe, 2*Fe + 1)]
        me_array = [me_array; Fe:-1:-Fe]
    end
    return hcat(Fe_array, me_array)
end

function fill_gDict(n2Fm_ats_g)
    inds = axes(n2Fm_ats_g, 1)
    return [repeat(inds, inner=length(inds)) repeat(inds, outer=length(inds))]
end

function fill_eDict(n2Fm_ats_e)
    inds = axes(n2Fm_ats_e, 1)
    return [repeat(inds, inner=length(inds)) repeat(inds, outer=length(inds))]
end

function mirror_around_end(array)
    return [array; 2*array[end] .- array[1:end-1] |> reverse]
end


""" Index - magnetic sublevel dictionaries """
function mF1(n2Fm_ats_g, g)
    return n2Fm_ats_g[g, 2]
end

function mF2(n2Fm_ats_e, e)
    return n2Fm_ats_e[e, 2]
end

function F1(n2Fm_ats_g, g)
    return n2Fm_ats_g[g, 1]
end

function F2(n2Fm_ats_e, e)
    return n2Fm_ats_e[e, 1]
end


"""
magn_field(par,B)

DESCRIPTION
Zeeman splitting H0B = H0 + HB
Computes magnetic sublevel energy splitting in magnetic field.

INPUT
par: atom parameters
B: magnetic field

OUTPUT
results_gr[1]: ground state eigenvalues
results_gr[2]: ground state eigenvectors
results_ex[1]: excited state eigenvalues
results_ex[2]: excited state eigenvectors

"""

function magn_field(par::param, laz::lazers, B, JIfmbasis_gr, Jzfmbasis_gr, Izfmbasis_gr, JIfmbasis_ex, Jzfmbasis_ex, Izfmbasis_ex)
    H = par.A_gr * JIfmbasis_gr + laz.muB_MHz * B * (par.gJ_gr * Jzfmbasis_gr + par.gI * Izfmbasis_gr)
    results_gr = eigenstates(dense((H + dagger(H)) / 2))

    H = par.A_ex * JIfmbasis_ex + laz.muB_MHz * B * (par.gJ_ex * Jzfmbasis_ex + par.gI * Izfmbasis_ex)
    results_ex = eigenstates(dense((H + dagger(H)) / 2))

    return results_gr[1], results_gr[2], results_ex[1], results_ex[2]
end


"""
DESCRIPTION
Transition matrix dipole elements - computed from Wigner 3j and 6j symbols

INPUT
F,m,J,I: atom quantum numbers
q: light polarization
g,e: magnetic sublevels corresponding to F, m quantum numbers. Called in rate equations

OUTPUT
Transition matrix dipole elements.

"""

function lin_three_j(f1, m1, f2, m2, q)
    if -m2 + m1 + q == 0
        value = (-1)^(f2 - m2) * wigner3j(f2, 1, f1, -m2, q, m1)
    else
        value = 0
    end
    return value
end

function lin_six_j(f1, f2, j1, j2, I)
    value = (-1)^(j2 + I + f1 + 1) * sqrt((2 * f1 + 1) * (2 * f2 + 1)) * wigner6j(j2, f2, I, f1, j1, 1)
    return value
end

function lin_three_j_star(f1, m1, f2, m2, q)
    value = (-1)^(f1 - m1) * wigner3j(f1, 1, f2, -m1, q, m2)
    return value
end

function lin_six_j_star(f1, f2, j1, j2, I)
    value = (-1)^(j1 + I + f2 + 1) * sqrt((2 * f1 + 1) * (2 * f2 + 1)) * wigner6j(j1, f1, I, f2, j2, 1)
    return value
end

function lin_dip(j1, j2, I, e, g, q, n2Fm_ats_g, n2Fm_ats_e)
    f1 = F1(n2Fm_ats_g, g)
    f2 = F2(n2Fm_ats_e, e)
    m1 = mF1(n2Fm_ats_g, g)
    m2 = mF2(n2Fm_ats_e, e)
    value = lin_six_j(f1, f2, j1, j2, I) * lin_three_j(f1, m1, f2, m2, q)
    return value
end

function lin_dip_star(j1, j2, I, e, g, q, n2Fm_ats_g, n2Fm_ats_e)
    f1 = F1(n2Fm_ats_g, g)
    f2 = F2(n2Fm_ats_e, e)
    m1 = mF1(n2Fm_ats_g, g)
    m2 = mF2(n2Fm_ats_e, e)
    value = lin_three_j_star(f1, m1, f2, m2, q) * lin_six_j_star(f1, f2, j1, j2, I)
    return value
end

function lmx_fillDipoleMatrix(par::param,e_vec,eigvects_gr,eigvects_ex,n2Fm_ats_g,n2Fm_ats_e)
    j1=par.J1
    j2=par.J2
    I=par.nucI
    dim_g=par.dim_g
    dim_e=par.dim_e

    dip_p=spzeros(Complex{Float64},dim_e,dim_g)
    for g in 1:1:dim_g
        for e in 1:1:dim_e
            dip_p[e,g]=e_vec[1]*lmx_dip(j1,j2,I,e,g,+1,eigvects_gr,eigvects_ex,n2Fm_ats_g,n2Fm_ats_e) 
        end
    end
    dip_0=spzeros(Complex{Float64},dim_e,dim_g)
    for g in 1:1:dim_g
        for e in 1:1:dim_e
            dip_0[e,g]=e_vec[2]*lmx_dip(j1,j2,I,e,g,0,eigvects_gr,eigvects_ex,n2Fm_ats_g,n2Fm_ats_e)
        end
    end
    dip_m=spzeros(Complex{Float64},dim_e,dim_g)
    for g in 1:1:dim_g
        for e in 1:1:dim_e
            dip_m[e,g]=e_vec[3]*lmx_dip(j1,j2,I,e,g,-1,eigvects_gr,eigvects_ex,n2Fm_ats_g,n2Fm_ats_e)
        end
    end
    dip = dip_p + dip_0 + dip_m
    return dip 
end

function lmx_fillDipoleMatrix_star(par::param,e_vec,eigvects_gr,eigvects_ex,n2Fm_ats_g,n2Fm_ats_e)
    j1=par.J1
    j2=par.J2
    I=par.nucI
    dim_g=par.dim_g
    dim_e=par.dim_e

    dip_p_star=spzeros(Complex{Float64},dim_g,dim_e) 
    for g in 1:1:dim_g
        for e in 1:1:dim_e
            dip_p_star[g,e]=conj(e_vec[3])*(-1)*lmx_dip_star(j1,j2,I,e,g,1,eigvects_gr,eigvects_ex,n2Fm_ats_g,n2Fm_ats_e) 
        end
    end
    dip_0_star=spzeros(Complex{Float64},dim_g,dim_e)
    for g in 1:1:dim_g
        for e in 1:1:dim_e
            dip_0_star[g,e]=conj(e_vec[2])*(-1)^0*lmx_dip_star(j1,j2,I,e,g,0,eigvects_gr,eigvects_ex,n2Fm_ats_g,n2Fm_ats_e)
        end
    end
    dip_m_star=spzeros(Complex{Float64},dim_g,dim_e)
    for g in 1:1:dim_g
        for e in 1:1:dim_e
            dip_m_star[g,e]=conj(e_vec[1])*(-1)*lmx_dip_star(j1,j2,I,e,g,-1,eigvects_gr,eigvects_ex,n2Fm_ats_g,n2Fm_ats_e)
        end
    end
    dip_star = dip_p_star + dip_0_star + dip_m_star
    return dip_star 
end

function lmx_dip(j1,j2,I,e,g,q,eigvects_gr,eigvects_ex,n2Fm_ats_g,n2Fm_ats_e)
    #f1=F1(g)
    #f2=F2(e)
    #m1=mF1(g)
    #m2=mF2(e)
    eps=1.0e-3
    f1=[]; f2=[]
    m1=[]; m2=[];
    coef_gr=[];coef_ex=[]
    index=0
    for element in eigvects_gr[g].data
        index += 1
        if norm(element) > eps
            f=F1(n2Fm_ats_g,index) #getF_gr(eigvect,par)
            m=mF1(n2Fm_ats_g,index) #getm_gr(eigvect,par)
            push!(f1,f)
            push!(m1,m)
            push!(coef_gr,norm(element))
        end
    end
    index=0
    for element in eigvects_ex[e].data
        index +=1
        if norm(element) > eps
            f=F2(n2Fm_ats_e,index) #getF_ex(eigvect,par)
            m=mF2(n2Fm_ats_e,index) #getm_ex(eigvect,par)
            push!(f2,f)
            push!(m2,m)
            push!(coef_ex,norm(element))
        end
    end
    value=0
    for i in eachindex(m1) #eachindex(f1)
        for j in eachindex(m2) #eachindex(f2)
            value += coef_gr[i]*coef_ex[j]*lin_six_j(f1[i],f2[j],j1,j2,I)*lin_three_j(f1[i],m1[i],f2[j],m2[j],q)
        end
    end
    #value=lin_six_j(f1,f2,j1,j2,I)*lin_three_j(f1,m1,f2,m2,q)
    return value
end

function lmx_dip_star(j1,j2,I,e,g,q,eigvects_gr,eigvects_ex,n2Fm_ats_g,n2Fm_ats_e)
    #f1=F1(g)
    #f2=F2(e)
    #m1=mF1(g)
    #m2=mF2(e)
    eps=1.0e-3
    f1=[]; f2=[]
    m1=[]; m2=[];
    coef_gr=[];coef_ex=[]
    index=0
    for element in eigvects_gr[g].data
        index += 1
        if norm(element) > eps          #eigvect is just a coefficient, not a ket
            f=F1(n2Fm_ats_g,index) #getF_gr(eigvect,par)     #need to get the original F,mF value corresponding to that position. 
            m=mF1(n2Fm_ats_g, index) #getm_gr(eigvect,par)     # maybe don't need operators. Already done in dict?
            push!(f1,f)
            push!(m1,m)
            push!(coef_gr,norm(element))
        end
    end
    index=0
    for element in eigvects_ex[e].data
        index += 1
        if norm(element) > eps
            f=F2(n2Fm_ats_e,index) #getF_ex(eigvect,par)
            m=mF2(n2Fm_ats_e,index) #getm_ex(eigvect,par)
            push!(f2,f)
            push!(m2,m)
            push!(coef_ex,norm(element))
        end
    end
    value=0
    for i in eachindex(m1) #eachindex(f1)
        for j in eachindex(m2) #eachindex(f2)
            value += coef_gr[i]*coef_ex[j]*lin_six_j_star(f1[i],f2[j],j1,j2,I)*lin_three_j_star(f1[i],m1[i],f2[j],m2[j],q)
        end
    end
    #value = lin_three_j_star(f1,m1,f2,m2,q)*lin_six_j_star(f1,f2,j1,j2,I)
    return value
end

function getF_gr(eigen_vector,par)
    F²_gr=par.F²fmbasis_gr
    F²=dagger(eigen_vector)*F²_gr*eigen_vector
    return round(-1+sqrt(1+4*F²)/2)
end
function getF_ex(eigen_vector,par)
    F²_ex=par.F²fmbasis_ex
    F²=dagger(eigen_vector)*F²_ex*eigen_vector
    return round(-1+sqrt(1+4*F²)/2)
end

function getm_gr(eigen_vector,par)
    m=dagger(eigen_vector)*par.Fzfmbasis_gr*eigen_vector
    return round(m)
end

function getm_ex(eigen_vector,par)
    m=dagger(eigen_vector)*par.Fzfmbasis_ex*eigen_vector
    return round(m)
end


""" 
Evals_g(Iz_gr,Jz_gr,F²fmbasis_gr,dim_g,eigvals_gr,eigvects_gr); Evals_e(Iz_ex,Jz_ex,F²fmbasis_ex,dim_e,eigvals_ex,eigvects_ex)

DESCRIPTION
Dictionary that connects calculated Zeeman splitting energy of a magnetic sublevel to F and m quantum numbers of energy levels.

INPUT
eigvals_gr: magnetic sublevel energies

OUTPUT
Dictionary for Fm quantum numbers of the level to the energy of the level.

"""

""" 
omega_g(par,Fm2E_g,gi,gj); omega_e(par,Fm2E_e,ei,ej)

DESCRIPTION
Computes energy difference between magnetic sublevels in ground state; in excited state.

INPUT
par: atom parameters
Fm2E: Fm to energy dictionary
gi,gj; ei,ej: index of magnetic sublevels         

OUTPUT
Energy difference between magnetic sublevels

"""

""" Energy difference between GROUND STATE levels """

function Evals_g(Iz_gr, Jz_gr, F²fmbasis_gr, dim_g, eigvals_gr, eigvects_gr)
    F²_gr = F²fmbasis_gr
    Fm2E_g = Dict{Tuple{Float64, Float64}, Float64}()


    for i = 1:dim_g

        Fz_gr = Iz_gr + Jz_gr

        mFg = dagger(eigvects_gr[i]) * Fz_gr * eigvects_gr[i]

        if round(real(mFg)) == -0.0
            mFg = 0.0
        end

        F² = dagger(eigvects_gr[i]) * F²_gr * eigvects_gr[i]
        #Fg = (-1 + sqrt(1 + 4 * F²)) / 2              ### FHG 2024-02-19
        Fg = -1+sqrt(1+4*F²)/2                         ### FHG 2024-02-19

        get!(Fm2E_g, (round(real(Fg)), round(real(mFg))), eigvals_gr[i])

    end

    return Fm2E_g

end

function omega_g(Fm2E_g, n2Fm_ats_g, gi, gj)
    return Fm2E_g[n2Fm_ats_g[gi, 1], n2Fm_ats_g[gi, 2]] - Fm2E_g[n2Fm_ats_g[gj, 1], n2Fm_ats_g[gj, 2]]
end

""" Energy difference between EXCITED STATE levels """
function Evals_e(Iz_ex, Jz_ex, F²fmbasis_ex, dim_e, eigvals_ex, eigvects_ex)
    F²_ex = F²fmbasis_ex
    Fm2E_e = Dict{Tuple{Float64, Float64}, Float64}()


    for i = 1:dim_e

        Fz_ex = Iz_ex + Jz_ex
        mFe = dagger(eigvects_ex[i]) * Fz_ex * eigvects_ex[i]
        if round(real(mFe)) == -0.0
            mFe = 0.0
        end
        F² = dagger(eigvects_ex[i]) * F²_ex * eigvects_ex[i]
        #Fe = (-1 + sqrt(1 + 4 * F²)) / 2           ### FHG 2024-02-19
        Fe = -1+sqrt(1+4*F²)/2                      ### FHG 2024-02-19   

        get!(Fm2E_e, (round(real(Fe)), round(real(mFe))), eigvals_ex[i])

    end

    return Fm2E_e

end

function omega_e(Fm2E_e, n2Fm_ats_e, ei, ej)
    return Fm2E_e[n2Fm_ats_e[ei, 1], n2Fm_ats_e[ei, 2]] - Fm2E_e[n2Fm_ats_e[ej, 1], n2Fm_ats_e[ej, 2]]
end

""" Energy difference between GROUND STATE and EXCITED STATE levels """
function omega_ge(par::param, g, e, Fm2E_g, Fm2E_e, n2Fm_ats_g, n2Fm_ats_e)
    par.Efs + Fm2E_e[n2Fm_ats_e[e, 1], n2Fm_ats_e[e, 2]] - Fm2E_g[n2Fm_ats_g[g, 1], n2Fm_ats_g[g, 2]]
end

"""
ksi_f(par::param,laz::lazers,g,e,Fm2E_g,Fm2E_e,detune) and ksi_cc_f(par::param,laz::lazers,e,g,Fm2E_g,Fm2E_e,detune)

DESCRIPTION
Creates Ξ matrix elements and fills them in a matrix.

INPUT
g: ground state level index, iterated through in rate equations
e: excited state level index, iterated through in rate equations
Fm2E_g: Fm to E ground state level dictionary 
Fm2E_e: Fm to E excited state level dictionary
detune: detuning (shift in frequency) due to Dopler effect

OUTPUT
Ξ matrix element

"""

#function ksi_f(par::param, laz::lazers, g, e, Fm2E_g, Fm2E_e, detune, n2Fm_ats_g, n2Fm_ats_e)
function ksi_f(par::param, laz::lazers, g, e, detune, eigvals_g, eigvals_e)
    #laz.Ωᵣ^2 / (((laz.Γ + laz.Δω) / 2) + laz.γ + 1.0im * (detune - omega_ge(par, g, e, Fm2E_g, Fm2E_e, n2Fm_ats_g, n2Fm_ats_e)))
    laz.Ωᵣ^2 / (((laz.Γ + laz.Δω) / 2) + laz.γ + 1.0im * (detune - (par.Efs + eigvals_e[e]-eigvals_g[g])))
end

#function ksi_cc_f(par::param, laz::lazers, e, g, Fm2E_g, Fm2E_e, detune, n2Fm_ats_g, n2Fm_ats_e)
function ksi_cc_f(par::param, laz::lazers, e, g, detune, eigvals_g, eigvals_e)
    #laz.Ωᵣ^2 / (((laz.Γ + laz.Δω) / 2) + laz.γ - 1.0im * (detune - omega_ge(par, g, e, Fm2E_g, Fm2E_e, n2Fm_ats_g, n2Fm_ats_e)))
    laz.Ωᵣ^2 / (((laz.Γ + laz.Δω) / 2) + laz.γ - 1.0im * (detune - (par.Efs + eigvals_e[e]-eigvals_g[g])))
end


"""
fill_ksi(par::param,laz::lazers,Fm2E_g,Fm2E_e,dim_g,dim_e,detune); fill_ksi_cc(par::param,laz::lazers,Fm2E_g,Fm2E_e,dim_e,dim_g,detune)

DESCRIPTION
Fills Ξ elements in a matrix.

INPUT
Fm2E_g: Fm to ground state energy dictionary 
Fm2E_e: Fm to excited state energy dictionary
dim_g: number of ground state levels
dim_e: number of excited state levels
detune: detuning (shift in frequency) due to Dopler effectv

OUTPUT
Ξ matrix; Ξ matrix complex conjugate

"""

#function fill_ksi(par::param, laz::lazers, Fm2E_g, Fm2E_e, dim_g, dim_e, detune, n2Fm_ats_g, n2Fm_ats_e)
function fill_ksi(par::param, laz::lazers, dim_g, dim_e, detune, eigvals_g, eigvals_e)
    ksi = spzeros(Complex{Float64}, dim_g, dim_e)
    for e = 1:dim_e
        for g = 1:dim_g
            #ksi[g, e] = ksi_f(par, laz, g, e, Fm2E_g, Fm2E_e, detune, n2Fm_ats_g, n2Fm_ats_e)
            ksi[g, e] = ksi_f(par, laz, g, e, detune, eigvals_g, eigvals_e)
        end
    end
    return ksi
end

#function fill_ksi_cc(par::param, laz::lazers, Fm2E_g, Fm2E_e, dim_e, dim_g, detune, n2Fm_ats_g, n2Fm_ats_e)
function fill_ksi_cc(par::param, laz::lazers, dim_e, dim_g, detune, eigvals_g, eigvals_e)    
    ksi_cc = spzeros(Complex{Float64}, dim_e, dim_g)
    for g = 1:dim_g
        for e = 1:dim_e
            #ksi_cc[e, g] = ksi_cc_f(par, laz, e, g, Fm2E_g, Fm2E_e, detune, n2Fm_ats_g, n2Fm_ats_e)
            ksi_cc[e, g] = ksi_cc_f(par, laz, e, g, detune, eigvals_g, eigvals_e)
        end
    end
    return ksi_cc
end


"""
GAMMA_f(Γ,g1,g2,e1,e2,j1,j2,I) 

DESCRIPTION
Calculates Γ matrix elements. Takes into account that two transitons should be excited with the same photon (mF2(e1)-mF1(g1) == mF2(e2)-mF1(g2)).

INPUT
Γ: spontaneous relaxation constant
g1,g2: ground state level index, iterated through in rate equations 
e1,e2: excited state level index, iterated through in rate equations
j1: quantum number of ground state total electron angular momentum
j2: quantum number of excited state total electron angular momentum
I: quantum number of total nuclear angular momentum

OUTPUT
Γ matrix element
"""

function GAMMA_f(Γ, g1, g2, e1, e2, j1, j2, I, n2Fm_ats_g, n2Fm_ats_e)
    sum = 0
    if mF2(n2Fm_ats_e, e1) - mF1(n2Fm_ats_g, g1) == mF2(n2Fm_ats_e, e2) - mF1(n2Fm_ats_g, g2)
        for q = -1:1
            sum += lin_dip(j1, j2, I, e1, g1, q, n2Fm_ats_g, n2Fm_ats_e) * lin_dip_star(j1, j2, I, e2, g2, -q, n2Fm_ats_g, n2Fm_ats_e) * (-1)^q
        end
    end
    sum *= (2 * j2 + 1) * Γ
    return sum
end


#function GAMMA_f_lmx(Γ,g1,g2,e1,e2,j1,j2,I,eigvects_gr, eigvects_ex,n2Fm_ats_g,n2Fm_ats_e) 
#    sum = 0
#    for i in 1:length(eigvects_gr[g1].data) 
#        for j in 1:length(eigvects_ex[e1].data)
#            for l in 1:length(eigvects_gr[g2].data)
#                for k in 1:length(eigvects_ex[e2].data)
#                    if mF2(n2Fm_ats_e,j)-mF1(n2Fm_ats_g,i) == mF2(n2Fm_ats_e,k)-mF1(n2Fm_ats_g,l)
##    if mF2(e1)-mF1(g1) == mF2(e2)-mF1(g2) 
#                       for q in -1:+1:+1
#                            #coeff=eig
#                            sum+=lmx_dip(j1,j2,I,e1,g1,q,eigvects_gr,eigvects_ex,n2Fm_ats_g,n2Fm_ats_e)*lmx_dip_star(j1,j2,I,e2,g2,-q,eigvects_gr,eigvects_ex,n2Fm_ats_g,n2Fm_ats_e)*(-1)^q
#                       end
#                    end
#                end
#            end
#        end    
#    end 
#    sum*=(2*j2+1)*Γ
#    return sum
#end

function GAMMA_f_lmx_2(Γ,g1,g2,e1,e2,j1,j2,I,eigvects_gr, eigvects_ex,n2Fm_ats_g,n2Fm_ats_e) 
    sum = 0
#    for i in 1:length(eigvects_gr[g1].data) 
#        for j in 1:length(eigvects_ex[e1].data)
#            for l in 1:length(eigvects_gr[g2].data)
#                for k in 1:length(eigvects_ex[e2].data)
                    #if mF2(n2Fm_ats_e,j)-mF1(n2Fm_ats_g,i) == mF2(n2Fm_ats_e,k)-mF1(n2Fm_ats_g,l)
                        if mF2(n2Fm_ats_e, e1) - mF1(n2Fm_ats_g, g1) == mF2(n2Fm_ats_e, e2) - mF1(n2Fm_ats_g, g2)
                        #    if mF2(e1)-mF1(g1) == mF2(e2)-mF1(g2) 
                       for q in -1:+1:+1
                            #coeff=eig
                            sum+=lmx_dip(j1,j2,I,e1,g1,q,eigvects_gr,eigvects_ex,n2Fm_ats_g,n2Fm_ats_e)*lmx_dip_star(j1,j2,I,e2,g2,-q,eigvects_gr,eigvects_ex,n2Fm_ats_g,n2Fm_ats_e)*(-1)^q
                       end
                    end
 #               end
 #           end
 #       end    
 #   end 
    sum*=(2*j2+1)*Γ
    return sum
end


"""
DESCRIPTION
Fills Γ matrix elements in a 4 dimensional matrix, where each element is described by 4 indexes of ground and excited state levels of two magnetic sublevel tranitions

INPUT
par: atom parameters
laz: laser parameters

OUTPUT
Γ matrix
"""
function fill_GAMMA(par::param, laz::lazers, n2Fm_ats_g, n2Fm_ats_e)
    GAMMA = zeros(Complex{Float64}, par.dim_g, par.dim_g, par.dim_e, par.dim_e)
    for e2 = 1:par.dim_e
        for e1 = 1:par.dim_e
            for g2 = 1:par.dim_g
                for g1 = 1:par.dim_g
                    GAMMA[g1, g2, e1, e2] = GAMMA_f(laz.Γ, g1, g2, e1, e2, par.J1, par.J2, par.nucI, n2Fm_ats_g, n2Fm_ats_e)
                end
            end
        end
    end
    return GAMMA
end

function fill_GAMMA_lmx(par::param,laz::lazers,eigvects_gr, eigvects_ex,n2Fm_ats_g,n2Fm_ats_e)
    GAMMA=zeros(Complex{Float64},par.dim_g,par.dim_g,par.dim_e,par.dim_e)
    for g1 in 1:1:par.dim_g
        for g2 in 1:1:par.dim_g
            for e1 in 1:1:par.dim_e 
                for e2 in 1:1:par.dim_e
                    GAMMA[g1,g2,e1,e2]=GAMMA_f_lmx_2(laz.Γ,g1,g2,e1,e2,par.J1,par.J2,par.nucI,eigvects_gr,eigvects_ex,n2Fm_ats_g,n2Fm_ats_e) 
                end
            end
        end
    end
    return GAMMA

end
"""
lin_fillDipoleMatrix(par::param,e_vec) and lin_fillDipoleMatrix_star(par::param,e_vec)

DESCRIPTION
Calculates dipole matrix elements using excited light vector components and fills them in a matrix. The dimensions of the matrix are the number of excited state levels x the number of ground state levels

INPUT
par: atom parameters
e_vec: excitation light E vector 

OUTPUT
Transition dipole matrix and its complex conjugate matrix.

"""

function lin_fillDipoleMatrix(par::param, e_vec, n2Fm_ats_g, n2Fm_ats_e)
    j1 = par.J1
    j2 = par.J2
    I = par.nucI
    dim_g = par.dim_g
    dim_e = par.dim_e

    dip_p = spzeros(Complex{Float64}, dim_e, dim_g)
    dip_0 = spzeros(Complex{Float64}, dim_e, dim_g)
    dip_m = spzeros(Complex{Float64}, dim_e, dim_g)
    for g = 1:dim_g
        for e = 1:dim_e
            dip_p[e, g] = e_vec[1] * lin_dip(j1, j2, I, e, g, +1, n2Fm_ats_g, n2Fm_ats_e)
            dip_0[e, g] = e_vec[2] * lin_dip(j1, j2, I, e, g, 0, n2Fm_ats_g, n2Fm_ats_e)
            dip_m[e, g] = e_vec[3] * lin_dip(j1, j2, I, e, g, -1, n2Fm_ats_g, n2Fm_ats_e)
        end
    end

    dip = dip_p + dip_0 + dip_m
    return dip
end

function lin_fillDipoleMatrix_star(par::param, e_vec, n2Fm_ats_g, n2Fm_ats_e)
    j1 = par.J1
    j2 = par.J2
    I = par.nucI
    dim_g = par.dim_g
    dim_e = par.dim_e

    dip_p_star = spzeros(Complex{Float64}, dim_g, dim_e)
    dip_0_star = spzeros(Complex{Float64}, dim_g, dim_e)
    dip_m_star = spzeros(Complex{Float64}, dim_g, dim_e)
    for e = 1:dim_e
        for g = 1:dim_g
            dip_p_star[g, e] = conj(e_vec[3]) * (-1) * lin_dip_star(j1, j2, I, e, g, 1, n2Fm_ats_g, n2Fm_ats_e)
            dip_0_star[g, e] = conj(e_vec[2]) * (-1)^0 * lin_dip_star(j1, j2, I, e, g, 0, n2Fm_ats_g, n2Fm_ats_e)
            dip_m_star[g, e] = conj(e_vec[1]) * (-1) * lin_dip_star(j1, j2, I, e, g, -1, n2Fm_ats_g, n2Fm_ats_e)
        end
    end

    dip_star = dip_p_star + dip_0_star + dip_m_star
    return dip_star
end


"""
ρDotgg(laz::lazers,matr,dim_g,dim_e,ksi,ksi_cc, dip, dip_star,Fm2E_g)

DESCRIPTION
Creates rate equation ground state coeficient matrix C:
dρ/dt=Cρ
Each row of the matrix represents one dρ(gigj)/dt equation, with each column representing coeffitient for each ρ element. This gives the matrix a dimension of (number of ground states)²x(number of ground states)². 

OUTPUT
Ground state coeficient matrix C.

"""

#function ρDotgg(γ, gDict, matr, dim_g, dim_e, ksi, ksi_cc, dip, dip_star, Fm2E_g, GAMMA, n2Fm_ats_g)
function ρDotgg(γ, gDict, matr, dim_g, dim_e, ksi, ksi_cc, dip, dip_star, eigvals_g, GAMMA)
    #for index = axes(gDict, 1)
    for gi = 1:dim_g  #gDict[index, 1]
      for  gj = 1:dim_g #gDict[index, 2]
        index=(gi-1)*dim_g+gj
        for ek = range(1, length=dim_e)
            for em = range(1, length=dim_e)
                if dip_star[gi, ek] * dip[em, gj] != 0
                    matr[index, dim_g*dim_g + (ek - 1)*dim_e + em] += (ksi[gi, em] + ksi_cc[ek, gj]) * dip_star[gi, ek] * dip[em, gj]
                end
            end

            for gm = range(1, length=dim_g)
                if dip_star[gi, ek] * dip[ek, gm] != 0
                    matr[index, (gm - 1)*dim_g + gj] += -ksi_cc[ek, gj] * dip_star[gi, ek] * dip[ek, gm]
                end

                if dip_star[gm, ek] * dip[ek, gj] != 0
                    matr[index, (gi - 1)*dim_g + gm] += -ksi[gi, ek] * dip_star[gm, ek] * dip[ek, gj]
                end
            end

            for el = range(1, length=dim_e)
                matr[index, dim_g*dim_g + (ek - 1)*dim_e + el] += GAMMA[gi, gj, ek, el]
            end
        end

        #matr[index, (gi - 1)*16 + gj] += -1.0im * omega_g(Fm2E_g, n2Fm_ats_g, gi, gj)
        matr[index, (gi - 1)*dim_g + gj] += -1.0im * (eigvals_g[gi]-eigvals_g[gj])
        matr[index, (gi - 1)*dim_g + gj] += -γ
      end
    end
end

"""
ρDotee(laz::lazers,matr,dim_g,dim_e,ksi,ksi_cc,dip, dip_star,Fm2E_e)

DESCRIPTION
Creates rate equation excitedd state coeficient matrix C:
dρ/dt=Cρ
Each row of the matrix represents one dρ(eiej)/dt equation, with each column representing coeffitient for each ρ element. This gives the matrix a dimension of (number of excited states)²x(number of excited states)². 

OUTPUT
Excited state coeficient matrix C.

"""

#function ρDotee(laz::lazers, eDict, gDict, matr, dim_g, dim_e, ksi, ksi_cc, dip, dip_star, Fm2E_e, n2Fm_ats_e)
function ρDotee(laz::lazers, eDict, gDict, matr, dim_g, dim_e, ksi, ksi_cc, dip, dip_star, eigvals_e)
    #for index = axes(eDict, 1)
    for ei = 1:dim_e #eDict[index, 1]
      for ej = 1:dim_e #eDict[index, 2]
        index=(ei-1)*dim_e+ej
        for gk = range(1, length=dim_g)
            for gm = range(1, length=dim_g)
                if dip_star[gm, ej] * dip[ei, gk] != 0
                    matr[dim_g*dim_g + index, (gk - 1)*dim_g + gm] += (ksi[gk, ej] + ksi_cc[ei, gm]) * dip_star[gm, ej] * dip[ei, gk]
                end
            end

            for em = range(1, length=dim_e)
                if dip_star[gk, em] * dip[ei, gk] != 0
                    matr[dim_g*dim_g + index, dim_g*dim_g + (em - 1)*dim_e + ej] += -ksi[gk, ej] * dip_star[gk, em] * dip[ei, gk]
                end

                if dip_star[gk, ej] * dip[em, gk] != 0
                    matr[dim_g*dim_g + index, dim_g*dim_g + (ei - 1)*dim_e + em] += -ksi_cc[ei, gk] * dip_star[gk, ej] * dip[em, gk]
                end
            end
        end

        #matr[size(gDict, 1) + index, size(gDict, 1) + (ei - 1)*16 + ej] += -1.0im * omega_e(Fm2E_e, n2Fm_ats_e, ei, ej)
        matr[dim_g*dim_g + index, dim_g*dim_g + (ei - 1)*dim_e + ej] += -1.0im * (eigvals_e[ei]-eigvals_e[ej])
        matr[dim_g*dim_g + index, dim_g*dim_g + (ei - 1)*dim_e + ej] += -(laz.Γ + laz.γ)
      end
    end
end


"""
nov_fillDipoleMatrix(j1,j2,I,dim_g,dim_e,e_vec) and nov_fillDipoleMatrix_star(j1,j2,I,dim_g,dim_e,e_vec).

DESCRIPTION
Calculates dipole matrix elements for observing using observation light vector components and fills them in a matrix. The dimensions of the matrix are the number of excited state levels x the number of ground state levels

INPUT
par: atom parameters
e_vec: observation light E vector 

OUTPUT
Transition dipole matrix and its complex conjugate matrix.

"""

#function nov_fillDipoleMatrix(par::param, e_vec, n2Fm_ats_g, n2Fm_ats_e)
function nov_fillDipoleMatrix(par::param, e_vec, eigvects_g, eigvects_e, n2Fm_ats_g, n2Fm_ats_e)
    j1 = par.J1
    j2 = par.J2
    I = par.nucI
    dim_g = par.dim_g
    dim_e = par.dim_e

    dip_p = spzeros(Complex{Float64}, dim_e, dim_g)
    dip_0 = spzeros(Complex{Float64}, dim_e, dim_g)
    dip_m = spzeros(Complex{Float64}, dim_e, dim_g)
    #for g = axes(n2Fm_ats_g, 1)
    #    for e = axes(n2Fm_ats_e, 1)
    #        dip_p[e, g] = e_vec[1] * lin_dip(j1, j2, I, e, g, +1, n2Fm_ats_g, n2Fm_ats_e)
    #        dip_0[e, g] = e_vec[2] * lin_dip(j1, j2, I, e, g, 0, n2Fm_ats_g, n2Fm_ats_e)
    #        dip_m[e, g] = e_vec[3] * lin_dip(j1, j2, I, e, g, -1, n2Fm_ats_g, n2Fm_ats_e)
    #    end
    #end
    for g in 1:1:dim_g  #g in keys(n2Fm_ats_g) 
        for e in 1:1:dim_e #e in keys(n2Fm_ats_e) 
            #dip_p[eDictInv[5,e],gDictInv[5,g]]=e_vec[1]*lin_dip(j1,j2,I,e,g,+1) 
            #dip_p[eDictInv[5,e],gDictInv[5,g]]=e_vec[1]*lmx_dip(j1,j2,I,e,g,+1,eigvects_gr,eigvects_ex,n2Fm_ats_g,n2Fm_ats_e) 
            #dip_0[eDictInv[5,e],gDictInv[5,g]]=e_vec[2]*lmx_dip(j1,j2,I,e,g,0,eigvects_g,eigvects_e,n2Fm_ats_g,n2Fm_ats_e)
            #dip_m[eDictInv[5,e],gDictInv[5,g]]=e_vec[3]*lmx_dip(j1,j2,I,e,g,-1,eigvects_g,eigvects_e,n2Fm_ats_g,n2Fm_ats_e)
            dip_p[e,g]=e_vec[1]*lmx_dip(j1,j2,I,e,g,+1,eigvects_g,eigvects_e,n2Fm_ats_g,n2Fm_ats_e) 
            dip_0[e,g]=e_vec[2]*lmx_dip(j1,j2,I,e,g,0,eigvects_g,eigvects_e,n2Fm_ats_g,n2Fm_ats_e)
            dip_m[e,g]=e_vec[3]*lmx_dip(j1,j2,I,e,g,-1,eigvects_g,eigvects_e,n2Fm_ats_g,n2Fm_ats_e)
        end
    end



    dip = dip_p + dip_0 + dip_m
    return dip
end

#function nov_fillDipoleMatrix_star(par::param, e_vec, n2Fm_ats_g, n2Fm_ats_e)
function nov_fillDipoleMatrix_star(par::param, e_vec, eigvects_g, eigvects_e, n2Fm_ats_g, n2Fm_ats_e)
    j1 = par.J1
    j2 = par.J2
    I = par.nucI
    dim_g = par.dim_g
    dim_e = par.dim_e

    dip_p_star = spzeros(Complex{Float64}, dim_g, dim_e)
    dip_0_star = spzeros(Complex{Float64}, dim_g, dim_e)
    dip_m_star = spzeros(Complex{Float64}, dim_g, dim_e)
    #for e = axes(n2Fm_ats_e, 1)
    #    for g = axes(n2Fm_ats_g, 1)
    #        dip_p_star[g, e] = conj(e_vec[3]) * (-1) * lin_dip_star(j1, j2, I, e, g, 1, n2Fm_ats_g, n2Fm_ats_e)
    #        dip_0_star[g, e] = conj(e_vec[2]) * (-1)^0 * lin_dip_star(j1, j2, I, e, g, 0, n2Fm_ats_g, n2Fm_ats_e)
    #        dip_m_star[g, e] = conj(e_vec[1]) * (-1) * lin_dip_star(j1, j2, I, e, g, -1, n2Fm_ats_g, n2Fm_ats_e)
    #    end
    #end
    for g in 1:1:dim_g  #g in keys(n2Fm_ats_g)
        for e in 1:1:dim_e #e in keys(n2Fm_ats_e)
            #dip_p_star[gDictInv[5,g],eDictInv[5,e]]=conj(e_vec[3])*(-1)*lin_dip_star(j1,j2,I,e,g,1)
            #dip_p_star[gDictInv[5,g],eDictInv[5,e]]=conj(e_vec[3])*(-1)*lmx_dip_star(j1,j2,I,e,g,1,eigvects_g,eigvects_e,n2Fm_ats_g,n2Fm_ats_e)
            #dip_0_star[gDictInv[5,g],eDictInv[5,e]]=conj(e_vec[2])*(-1)^0*lmx_dip_star(j1,j2,I,e,g,0,eigvects_g,eigvects_e,n2Fm_ats_g,n2Fm_ats_e)
            #dip_m_star[gDictInv[5,g],eDictInv[5,e]]=conj(e_vec[1])*(-1)*lmx_dip_star(j1,j2,I,e,g,-1,eigvects_g,eigvects_e,n2Fm_ats_g,n2Fm_ats_e)
            dip_p_star[g,e]=conj(e_vec[3])*(-1)*lmx_dip_star(j1,j2,I,e,g,1,eigvects_g,eigvects_e,n2Fm_ats_g,n2Fm_ats_e)
            dip_0_star[g,e]=conj(e_vec[2])*(-1)^0*lmx_dip_star(j1,j2,I,e,g,0,eigvects_g,eigvects_e,n2Fm_ats_g,n2Fm_ats_e)
            dip_m_star[g,e]=conj(e_vec[1])*(-1)*lmx_dip_star(j1,j2,I,e,g,-1,eigvects_g,eigvects_e,n2Fm_ats_g,n2Fm_ats_e)
        end
    end


    dip_star = dip_p_star + dip_0_star + dip_m_star
    return dip_star
end






#dstep_length=dscan/2

#function signals(B₀, par, laz, gDict, eDict, n2Fm_ats_g, n2Fm_ats_e, dscan=dscan, dstep_length=dstep_length, e_vec_n=e_vec_n, e_vect_z=e_vect_z, sigma=sigma)
function signals(B₀, par, laz, evecs, Doppler_steps)
    # polarization vectors
    e_vec_i=evecs[1]
    e_vec_n=evecs[2]
    e_vect_z=evecs[3]

    # Dictionaries
    n2Fm_ats_g = fill_n2Fm_ats_g(par)
    n2Fm_ats_e = fill_n2Fm_ats_e(par)

    gDict = fill_gDict(n2Fm_ats_g)
    eDict = fill_eDict(n2Fm_ats_e)
    #induced relaxation
    λ = zeros(Complex{Float64}, par.dim_g * par.dim_g + par.dim_e * par.dim_e)
    for index = axes(gDict, 1)
        gi = gDict[index, 1]
        gj = gDict[index, 2]
        if (gi == gj)
            λ[index] = laz.γ / par.dim_g
        end
    end

    """Dopler shift"""
    step_lit = 2
    sigma = sqrt(laz.kB * laz.T / laz.masa) * laz.ω_svitr / laz.c
    dscan = 2 * sigma * step_lit
    dstep_nr = Doppler_steps * step_lit #number of steps
    dstep_length = dscan / dstep_nr

    println("Processing magnetic field B₀=", B₀, " G."      )
    eigvals_gr, eigvects_gr, eigvals_ex, eigvects_ex = magn_field(par, laz, B₀,
                                                                  par.JIfmbasis_gr, par.Jzfmbasis_gr, par.Izfmbasis_gr,
                                                                  par.JIfmbasis_ex, par.Jzfmbasis_ex, par.Izfmbasis_ex)

    
    #local dip_star = lmx_fillDipoleMatrix_star(par,e_vec_i,eigvects_gr,eigvects_ex,n2Fm_ats_g,n2Fm_ats_e) #izveido dipolu matricu
    #local dip=lmx_fillDipoleMatrix(par,e_vec_i,eigvects_gr,eigvects_ex,n2Fm_ats_g,n2Fm_ats_e)
    
    local dip_star = lin_fillDipoleMatrix_star(par, e_vec_i, n2Fm_ats_g, n2Fm_ats_e) #izveido dipolu matricu
    local dip = lin_fillDipoleMatrix(par, e_vec_i, n2Fm_ats_g, n2Fm_ats_e)

    #Fm2E_g = Evals_g(par.Izfmbasis_gr, par.Jzfmbasis_gr, par.F²fmbasis_gr, par.dim_g, eigvals_gr, eigvects_gr)
    #Fm2E_e = Evals_e(par.Izfmbasis_ex, par.Jzfmbasis_ex, par.F²fmbasis_ex, par.dim_e, eigvals_ex, eigvects_ex)

    GAMMA = fill_GAMMA_lmx(par,laz,eigvects_gr, eigvects_ex, n2Fm_ats_g,n2Fm_ats_e)

    A_Doplera = 0
    I_Doplera = 0

    local rho_e = spzeros(Complex{Float64}, par.dim_e * par.dim_e)
    local rho_g = spzeros(Complex{Float64}, par.dim_g * par.dim_g)   ### FHG 2024-02-19

    for dshift = range(-dscan, step=dstep_length, stop=dscan)
        detune = laz.ω_svitr - dshift

        matr = spzeros(Complex{Float64}, par.dim_g * par.dim_g + par.dim_e * par.dim_e, par.dim_g * par.dim_g + par.dim_e * par.dim_e)
        #ksi = fill_ksi(par, laz, Fm2E_g, Fm2E_e, par.dim_g, par.dim_e, detune, n2Fm_ats_g, n2Fm_ats_e)
        ksi = fill_ksi(par, laz, par.dim_g, par.dim_e, detune, eigvals_gr, eigvals_ex)
        #ksi_cc = fill_ksi_cc(par, laz, Fm2E_g, Fm2E_e, par.dim_e, par.dim_g, detune, n2Fm_ats_g, n2Fm_ats_e)
        ksi_cc = fill_ksi_cc(par, laz, par.dim_e, par.dim_g, detune, eigvals_gr, eigvals_ex)

        #ρDotgg(laz.γ, gDict, matr, par.dim_g, par.dim_e, ksi, ksi_cc, dip, dip_star, Fm2E_g, GAMMA, n2Fm_ats_g)
        ρDotgg(laz.γ, gDict, matr, par.dim_g, par.dim_e, ksi, ksi_cc, dip, dip_star, eigvals_gr, GAMMA)
        #ρDotee(laz, eDict, gDict, matr, par.dim_g, par.dim_e, ksi, ksi_cc, dip, dip_star, Fm2E_e, n2Fm_ats_e)
        ρDotee(laz, eDict, gDict, matr, par.dim_g, par.dim_e, ksi, ksi_cc, dip, dip_star, eigvals_ex)

        lu_matr_sparse = lu(matr)

        eigen_vectors = lu_matr_sparse \ (-λ)

        @views rho_e = permutedims(reshape(eigen_vectors[par.dim_g*par.dim_g+1:(par.dim_g*par.dim_g)+par.dim_e*par.dim_e], par.dim_e, par.dim_e))
        rho_e_dimge = spzeros(Complex{Float64}, par.dim_g + par.dim_e, par.dim_g + par.dim_e)
        @views rho_e_dimge[(par.dim_g+1):(par.dim_g+par.dim_e), (par.dim_g+1):(par.dim_g+par.dim_e)] = rho_e

        d_star_nov_dimge = spzeros(Complex{Float64}, par.dim_g + par.dim_e, par.dim_g + par.dim_e)
        #d_star_nov = nov_fillDipoleMatrix_star(par, e_vec_n, n2Fm_ats_g, n2Fm_ats_e)
        d_star_nov = nov_fillDipoleMatrix_star(par, e_vec_n, eigvects_gr, eigvects_ex, n2Fm_ats_g, n2Fm_ats_e)
        @views d_star_nov_dimge[1:par.dim_g, (par.dim_g+1):(par.dim_g+par.dim_e)] = d_star_nov

        d_nov_dimge = spzeros(Complex{Float64}, par.dim_g + par.dim_e, par.dim_g + par.dim_e)
        #d_nov = nov_fillDipoleMatrix(par, e_vec_n, n2Fm_ats_g, n2Fm_ats_e)
        d_nov = nov_fillDipoleMatrix(par, e_vec_n, eigvects_gr, eigvects_ex,n2Fm_ats_g, n2Fm_ats_e)
        @views d_nov_dimge[(par.dim_g+1):(par.dim_g+par.dim_e), 1:par.dim_g] = d_nov

        I = d_star_nov_dimge * rho_e_dimge * d_nov_dimge
        I_Doplera += tr(I) * dstep_length * (exp(-dshift^2 / (2 * sigma^2))) / (sqrt(2 * π) * sigma)

        @views rho_g = permutedims(reshape(eigen_vectors[1:par.dim_g*par.dim_g], par.dim_g, par.dim_g))

        GammaR2 = ((laz.Γ + laz.Δω + laz.γ) / 2)^2
        #d_star_z = nov_fillDipoleMatrix_star(par, e_vect_z, n2Fm_ats_g, n2Fm_ats_e)
        d_star_z = nov_fillDipoleMatrix_star(par, e_vect_z, eigvects_gr, eigvects_ex, n2Fm_ats_g, n2Fm_ats_e)
        #d_z = nov_fillDipoleMatrix(par, e_vect_z, n2Fm_ats_g, n2Fm_ats_e)
        d_z = nov_fillDipoleMatrix(par, e_vect_z, eigvects_gr, eigvects_ex, n2Fm_ats_g, n2Fm_ats_e)

        A = 0
        for gj = axes(n2Fm_ats_g, 1)
            for gk = axes(n2Fm_ats_g, 1)
                for ei = axes(n2Fm_ats_e, 1)
                    A += (d_z[ei, gj] * rho_g[gj, gk] * d_star_z[gk, ei]) /
                         #(GammaR2 + (detune - omega_ge(par, gj, ei, Fm2E_g, Fm2E_e, n2Fm_ats_g, n2Fm_ats_e))^2)
                         (GammaR2 + (detune - (par.Efs + eigvals_ex[ei] - eigvals_gr[gj]))^2)
                end
            end
        end
        A_Doplera += A * dstep_length * (exp(-(dshift^2) / (2 * sigma^2))) / (sqrt(2 * π) * sigma)
    end
    
    #return (I_Doplera, A_Doplera) # [I_Doplera, A_Doplera]
    return (I_Doplera, A_Doplera, rho_g, rho_e) # [I_Doplera, A_Doplera]   ### FHG 2024-02-19 added return for rho_g, rho_e
end

#signals_for_pmap(B₀) = signals(B₀, par, laz, gDict, eDict, n2Fm_ats_g, n2Fm_ats_e)




end
