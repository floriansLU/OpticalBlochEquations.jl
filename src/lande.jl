function LandeFactorJ(J,L,S)
    return 1+(J*(J+1)+S*(S+1)-L*(L+1))/(2*J*(J+1))
end



function LandeFactorF(F,I,J, gJ) 
    return gJ*(F*(F+1)-I*(I+1)+J*(J+1))/(2*F*(F+1))
end