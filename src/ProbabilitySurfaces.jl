"""Probability Surfaces"""

using FastTransforms
using QuantumOptics
using WignerSymbols

# to get spherical harmonic Ylm, use package FastTransforms:
# Ylm(θ,ϕ)=FastTransforms.sphevaluate(θ,ϕ,l,m).
function probSurface(J,ρ,θ,ϕ)
    norm=sqrt(4*pi/(2J*+1))
    sum=0
    for k in 0:1:2*J
        for q in -k:1:k
            Tqk=reshape(dagger(T(q,k,J)).data,Int64(2*J+1),Int64(2*J+1))
            ρqk=tr(ρ*Tqk)
            Ykq=sphevaluate(θ,ϕ,k,q)
            sum += clebschgordan(J,J,k,0,J,J)*ρqk*Ykq
        end
    end
    return sum
end

function T(q,k,J)
    Jbasis = SpinBasis(J)
    operator = DenseOperator(Jbasis,Jbasis,zeros(Int(2*J+1),Int(2*J+1)))
    lastState=spindown(Jbasis)
    mDict=Dict()
    for m in -J:1:J
        mDict[m] = lastState/norm(lastState)
        lastState=sigmap(Jbasis)*lastState
    end
    flag=0
    for m in -J:1:J
        if (q-m >= -J) && (q-m <= J) && (-q >= -k && -q) <= k 
            if flag==0 
                operator = (-1)^(J-m)*(2*k+1)^(1/2)*wigner3j(J,J,k,m,q-m,-q)*mDict[m] ⊗ mDict[m-q]
                flag=1
            else
                operator += (-1)^(J-m)*(2*k+1)^(1/2)*wigner3j(J,J,k,m,q-m,-q)*mDict[m] ⊗ mDict[m-q]
            end
        end
    end
    return operator
end

function TT(q,k,J)
    Jbasis = SpinBasis(J)
    operator = DenseOperator(Jbasis,Jbasis,zeros(Int(2*J+1),Int(2*J+1)))
    lastState=spindown(Jbasis)
    mDict=Dict()
    for m in -J:1:J
        mDict[m] = lastState/norm(lastState)
        lastState=sigmap(Jbasis)*lastState
    end
    flag=0
    for m1 in -J:1:J 
        for m2= -J:1:J
            if q==m1-m2 
                if flag==0 && clebschgordan(J,m1,J,-m2,k,q) != 0
                    operator = (-1)^(J-m2)*clebschgordan(J,m1,J,-m2,k,q)*mDict[m1] ⊗ mDict[m2]
                    flag=1
                else
                    operator += (-1)^(J-m2)*clebschgordan(J,m1,J,-m2,k,q)*mDict[m1] ⊗ mDict[m2]
                end
            end
        end
    end
    return operator
end

function Xcoord(J,ρ,θ,ϕ)
    r=real(probSurface(J,ρ,θ,ϕ))
    return r*sin(θ)*sin(ϕ)
end

function Ycoord(J,ρ,θ,ϕ)
    r=real(probSurface(J,ρ,θ,ϕ))
    return r*sin(θ)*cos(ϕ)
end

function Zcoord(J,ρ,θ,ϕ)
    r=real(probSurface(J,ρ,θ,ϕ))
    return r*cos(θ)
end

function plotProbSurf(J,ρ,θrange=0:0.1:pi,ϕrange=0:0.1:2*pi)
    xs=[Xcoord(J,ρ,θ,ϕ) for θ in θrange, ϕ in ϕrange ]
    ys=[Ycoord(J,ρ,θ,ϕ) for θ in θrange, ϕ in ϕrange ]
    zs=[Zcoord(J,ρ,θ,ϕ) for θ in θrange, ϕ in ϕrange ]
    surface(xs,ys,zs)
end