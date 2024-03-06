using QuantumOptics
using WignerSymbols
using Parameters

#include("../src/AtomicTransitions.jl")
include("../src/lande.jl")

"""
dim_g: number of ground state magnetic sublevels
dim_e: number of excited state magnetic sublevels
gJ_gr: ground state fine structure Lande factor
gJ_ex: excited state fine structure Lande factor
mJmIbasis_gr: mJmI composite basis
CG: Clebschgordan coeffitient matrix
"""
struct param

    J1::Float64
    J2::Float64
    nucI::Float64
    exFmax::Float64
    exFmin::Float64
    grFmax::Float64
    grFmin::Float64
    S1::Float64
    S2::Float64
    L1::Float64
    L2::Float64

    gI::Float64

    A_gr::Float64
    A_ex::Float64
    B_ex::Float64
    Efs::Float64


    dim_g::Int64
    dim_e::Int64
    gJ_gr::Float64
    gJ_ex::Float64

    Jzfmbasis_gr::Operator
    JIfmbasis_gr::Operator
    Jzfmbasis_ex::Operator
    JIfmbasis_ex::Operator
    mJmIbasis_ex::CompositeBasis
    Izfmbasis_gr::Operator
    Izfmbasis_ex::Operator

    Fzfmbasis_gr::Operator
    Fzfmbasis_ex::Operator
    F²fmbasis_gr::Operator
    F²fmbasis_ex::Operator
    Iz_gr::Operator   #### FHG 2024-02-19
    Jz_gr::Operator   ### FHG 2024-02-19
    Iz_ex::Operator
    Jz_ex::Operator
    mJmIbasisF²_ex::Operator

    function param(cezijsD)

        J1 = cezijsD().J1
        J2 = cezijsD().J2
        nucI = cezijsD().nucI
        exFmax = cezijsD().exFmax
        exFmin = cezijsD().exFmin
        grFmax = cezijsD().grFmax
        grFmin = cezijsD().grFmin
        S1 = cezijsD().S1
        S2 = cezijsD().S2
        L1 = cezijsD().L1
        L2 = cezijsD().L2

        gI = cezijsD().gI

        A_gr = cezijsD().A_gr
        A_ex = cezijsD().A_ex
        B_ex = cezijsD().B_ex
        Efs = cezijsD().Efs


        dim_g = Int64((2 * J1 + 1) * (2 * nucI + 1))
        dim_e = Int64((2 * J2 + 1) * (2 * nucI + 1))

        gJ_gr = LandeFactorJ(J1, L1, S1)
        gJ_ex = LandeFactorJ(J2, L2, S2)


        """ Ground state """
        mJbasis_gr = SpinBasis(J1)
        mIbasis = SpinBasis(nucI)
        mJmIbasis_gr = CompositeBasis(mJbasis_gr, mIbasis)


        mIbasisSz = 0.5 * sigmaz(mIbasis)
        mJbasisSz_gr = 0.5 * sigmaz(mJbasis_gr)
        mIbasisS₊ = sigmap(mIbasis)
        mIbasisS₋ = sigmam(mIbasis)
        mJbasisS₊_gr = sigmap(mJbasis_gr)
        mJbasisS₋_gr = sigmam(mJbasis_gr)
        mIbasisIdentity = identityoperator(mIbasis)
        mJbasisIdentity_gr = identityoperator(mJbasis_gr)



        Iz_gr = tensor(mJbasisIdentity_gr, mIbasisSz)
        Jz_gr = tensor(mJbasisSz_gr, mIbasisIdentity)
        I₊_gr = tensor(mJbasisIdentity_gr, mIbasisS₊)
        I₋_gr = tensor(mJbasisIdentity_gr, mIbasisS₋)
        J₊_gr = tensor(mJbasisS₊_gr, mIbasisIdentity)
        J₋_gr = tensor(mJbasisS₋_gr, mIbasisIdentity)
        JI_gr = Jz_gr * Iz_gr + 0.5 * (J₊_gr * I₋_gr + J₋_gr * I₊_gr)

        mJmIbasisI²_gr = nucI * (nucI + 1) * identityoperator(mJmIbasis_gr)
        mJmIbasisJ²_gr = J1 * (J1 + 1) * identityoperator(mJmIbasis_gr)
        mJmIbasisF²_gr = mJmIbasisI²_gr + mJmIbasisJ²_gr + 2 * Iz_gr * Jz_gr + I₊_gr * J₋_gr + I₋_gr * J₊_gr

        j1 = nucI  # I
        j2 = J1  # J
        row = 1
        col = 1
        CG_gr = zeros(Float64, dim_g, dim_g)
        for j3 = grFmax:-1:grFmin  #F
            for m3 = j3:-1:-j3  #mF
                for m1 = j1:-1:-j1  #mI 
                    for m2 = j2:-1:-j2  #mJ
                        CG_gr[row, col] = clebschgordan(j1, m1, j2, m2, j3, m3)
                        col += 1
                    end

                end
                row += 1
                col = 1
            end
        end

        CG_gr = Float64.(CG_gr)


        TransformMatrix_gr = DenseOperator(mJmIbasis_gr, CG_gr)
        TransformMatrixDagger_gr = dagger(TransformMatrix_gr)


        Jzfmbasis_gr = TransformMatrix_gr * (Jz_gr) * TransformMatrixDagger_gr
        Izfmbasis_gr = TransformMatrix_gr * (Iz_gr) * TransformMatrixDagger_gr
        JIfmbasis_gr = TransformMatrix_gr * (JI_gr) * TransformMatrixDagger_gr
        Fzfmbasis_gr = Jzfmbasis_gr + Izfmbasis_gr

        F²fmbasis_gr = TransformMatrix_gr * (mJmIbasisF²_gr) * TransformMatrixDagger_gr

        """ Excited state """

        mJbasis_ex = SpinBasis(J2)
        mIbasis_ex = SpinBasis(nucI)
        mJmIbasis_ex = CompositeBasis(mJbasis_ex, mIbasis_ex)

        mIbasisSz_ex = 0.5 * sigmaz(mIbasis_ex)
        mJbasisSz_ex = 0.5 * sigmaz(mJbasis_ex)
        mIbasisS₊_ex = sigmap(mIbasis_ex)
        mIbasisS₋_ex = sigmam(mIbasis_ex)
        mJbasisS₊_ex = sigmap(mJbasis_ex)
        mJbasisS₋_ex = sigmam(mJbasis_ex)
        mIbasisIdentity_ex = identityoperator(mIbasis_ex)
        mJbasisIdentity_ex = identityoperator(mJbasis_ex)


        Iz_ex = tensor(mJbasisIdentity_ex, mIbasisSz_ex)
        Jz_ex = tensor(mJbasisSz_ex, mIbasisIdentity_ex)
        I₊_ex = tensor(mJbasisIdentity_ex, mIbasisS₊_ex)
        I₋_ex = tensor(mJbasisIdentity_ex, mIbasisS₋_ex)
        J₊_ex = tensor(mJbasisS₊_ex, mIbasisIdentity_ex)
        J₋_ex = tensor(mJbasisS₋_ex, mIbasisIdentity_ex)
        JI_ex = Jz_ex * Iz_ex + 0.5 * (J₊_ex * I₋_ex + J₋_ex * I₊_ex)


        mJmIbasisI²_ex = nucI * (nucI + 1) * identityoperator(mJmIbasis_ex)
        mJmIbasisJ²_ex = J2 * (J2 + 1) * identityoperator(mJmIbasis_ex)
        mJmIbasisF²_ex = mJmIbasisI²_ex + mJmIbasisJ²_ex + 2 * Iz_ex * Jz_ex + I₊_ex * J₋_ex + I₋_ex * J₊_ex

        j1 = nucI  # I
        j2 = J2  # J
        row = 1
        col = 1
        CG_ex = zeros(Float64, dim_e, dim_e)
        for j3 = exFmax:-1:exFmin  #F
            for m3 = j3:-1:-j3  #mF
                for m1 = j1:-1:-j1  #mI
                    for m2 = j2:-1:-j2  #mJ
                        CG_ex[row, col] = clebschgordan(j1, m1, j2, m2, j3, m3)
                        col += 1
                    end
                end
                row += 1
                col = 1
            end
        end
        CG_ex


        CG_ex = Float64.(CG_ex)


        TransformMatrix_ex = DenseOperator(mJmIbasis_ex, CG_ex)
        TransformMatrixDagger_ex = dagger(TransformMatrix_ex)


        Jzfmbasis_ex = TransformMatrix_ex * (Jz_ex) * TransformMatrixDagger_ex
        Izfmbasis_ex = TransformMatrix_ex * (Iz_ex) * TransformMatrixDagger_ex
        JIfmbasis_ex = TransformMatrix_ex * (JI_ex) * TransformMatrixDagger_ex
        Fzfmbasis_ex = Jzfmbasis_ex + Izfmbasis_ex


        F²fmbasis_ex = TransformMatrix_ex * (mJmIbasisF²_ex) * TransformMatrixDagger_ex
        #new(J1, J2, nucI, exFmax, exFmin, grFmax, grFmin, S1, S2, L1, L2, gI, A_gr, A_ex, B_ex, Efs, dim_g, dim_e, gJ_gr, gJ_ex, Jzfmbasis_gr, JIfmbasis_gr, Jzfmbasis_ex, JIfmbasis_ex, mJmIbasis_ex, Izfmbasis_gr, Izfmbasis_ex, Fzfmbasis_gr, Fzfmbasis_ex, F²fmbasis_gr, F²fmbasis_ex, Iz_ex, Jz_ex, mJmIbasisF²_ex)  ### FHG 2024-02-19
        new(J1, J2, nucI, exFmax, exFmin, grFmax, grFmin, S1, S2, L1, L2, gI, A_gr, A_ex, B_ex, Efs, dim_g, dim_e, gJ_gr, gJ_ex, Jzfmbasis_gr, JIfmbasis_gr, Jzfmbasis_ex, JIfmbasis_ex, mJmIbasis_ex, Izfmbasis_gr, Izfmbasis_ex, Fzfmbasis_gr, Fzfmbasis_ex, F²fmbasis_gr, F²fmbasis_ex, Iz_gr, Jz_gr, Iz_ex, Jz_ex, mJmIbasisF²_ex)  ### FHG 2024-02-19

    end

end

"""
Ωᵣ: Rabi frequency
Γ: spontaneous relaxation constant
γ: induced relaxation
ω_svitr: central laser frequency
Δω: laser linewidth
kB: Bolcman constant
T: Temperature
masa: mass of one atom of the vapor
c: speed of light
muB_MHz: Bohr magneton

"""
@with_kw struct lazers
    Ωᵣ::ComplexF64 = 1.0 + 0.0im
    Γ::ComplexF64 = 4.575 + 0.0im
    γ::ComplexF64 = 0.019 + 0.0im
    ω_svitr::Float64 = 335120562.8#D1 3->3 transition
    Δω::Float64 = 4 #2 #MHz  #FHG 2024-01-31 
    kB::Float64 = 1.3806504 * 10^-23
    T::Float64 = 298
    masa::Float64 = 2.206211106 * 10^-25
    c::Float64 = 2.99792458 * 10^8
    muB_MHz::Float64 = 1.3996 #MHz
end

