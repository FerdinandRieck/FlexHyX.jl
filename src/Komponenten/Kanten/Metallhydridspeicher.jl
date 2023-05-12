Base.@kwdef mutable struct y_mMH
    m::Number = 0
    e::Number = 0
end

Base.@kwdef mutable struct mMH_kante <: Gas_Kante
    #-- default Parameter
    Param::GPMH_Param

    #-- Gasknoten links und rechts
    KL::Gas_Knoten
    KR::GPMH_Knoten

    #-- Zustandsvariablen
    y = y_mMH()

    #-- M-Matrix
    M::Array{Int} = [0; 0] 

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mMH_kante,t)
        #-- Parameter
        (; KAepsdivmuL, delta_rho,V_s,V_g,cv) = kante.Param
        #--
    
        #-- Zustandsvariablen 
        m = kante.y.m
        e = kante.y.e
        #--
    
        #-- Knoten links und rechts
        (; KL,KR) = kante
        PL = KL.y.P
        PR = KR.y.P
        TL = KL.y.T
        TR = KR.y.T
        Θ = KR.y.Θ
        M_H2 = KR.y.M_H2
        #--

        rho_s = delta_rho*Θ
        rho_g = (M_H2 - rho_s*V_s)/V_g
        delta_P = PL - PR;
        md = KAepsdivmuL*rho_g*delta_P;  #--- Massenfluss

        dy[k] = m - md
        dy[k+1] = e - cv*m*ifxaorb(m,TL,TR)
end
