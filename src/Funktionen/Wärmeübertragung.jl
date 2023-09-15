#-- dimensionslosen Kenngrößen
function Nusselt(Re,Pr,D,L)
    if typeof(Re) != Symbolics.Num
        if  Re < 2300  #-- laminare Strömung
            Nu = (49.37+(1.615*(Re*Pr*D/L)^(1/3)-0.7)^3)^(1/3)
        elseif Re < 1.0e4 #-- Übergangsbereich
            Nu_lam = (49.37+(1.615*(Re*Pr*D/L)^(1/3)-0.7)^3)^(1/3)
            ξ = (1.82*log10(Re)-1.64)^(-2)
            Nu_turb = ξ/8*(Re-1000)*Pr/(1+12.7*(ξ/8)^(1/2)*(Pr^(2/3)-1))*(1+(D/L)^(2/3))
            gamma = (Re - 2300)/(1.0e4 - 2300)
            Nu = (1-gamma)*Nu_lam + gamma*Nu_turb
        else #-- turbulente Strömung
            ξ = (1.82*log10(Re)-1.64)^(-2)
            Nu = ξ/8*(Re-1000)*Pr/(1+12.7*(ξ/8)^(1/2)*(Pr^(2/3)-1))*(1+(D/L)^(2/3))
        end
    else
        Nu_lam = (49.37+(1.615*(Re*Pr*D/L)^(1/3)-0.7)^3)^(1/3)
        ξ = (1.82*log10(Re)-1.64)^(-2)
        Nu_turb = ξ/8*(Re-1000)*Pr/(1+12.7*(ξ/8)^(1/2)*(Pr^(2/3)-1))*(1+(D/L)^(2/3))
        gamma = (Re - 2300)/(1.0e4 - 2300)
        Nu = (1-gamma)*Nu_lam + gamma*Nu_turb
    end
    return Nu
end

function Reynolds(m,D,A,mu)
    Re = abs(m)*D/(A*mu); Re = max(Re,1.0e-6)
    return Re
end

function Prandtl(cp,mu,lam)
    Pr = cp*mu/lam
    return Pr
end
#------------------------

function Wärmeübergang(Nu,lam,D)
    alpha = Nu*lam/D
    return alpha
end

function Wärmeduchgang(alpha_i,alpha_a,Di,Da,lamRohr)
    kA = (1/alpha_i*Da/Di + Da/2*(1/lamRohr*log(Da/Di)) + 1/alpha_a)^-1
    return kA
end

function Wärmeduchgang(alpha_i,Di,Dm,Da,lamRohr1,lamRohr2)
    kA = (1/alpha_i*Da/Di + Da/2*(1/lamRohr1*log(Dm/Di) + 1/lamRohr2*log(Da/Dm)))^-1
    return kA
end

function Wärmeduchgang(alpha_i,alpha_a,Di,Dm,Da,lamRohr1,lamRohr2)
    kA = (1/alpha_i*Da/Di + Da/2*(1/lamRohr1*log(Dm/Di) + 1/lamRohr2*log(Da/Dm)) + 1/alpha_a)^-1
    return kA
end

#-- Berechnung Wärmedurchgangskoeffizienten Rohr
function Wärmeduchgang_Rohr_aussen(kante)
    m = kante.y.m
    Di = kante.Param.Di
    Dm = kante.Param.Dm
    Da = kante.Param.Da
    A = kante.Param.A
    L = kante.Param.L
    mu = kante.Param.mu
    cv_H2O = kante.Param.cv_H2O
    lamW = kante.Param.lamW
    lamRohr1 = kante.Param.lamRohr1
    lamRohr2 = kante.Param.lamRohr2
    Re = Reynolds(m,Di,A,mu)
    Pr = Prandtl(cv_H2O,mu,lamW)
    Nu = Nusselt(Re,Pr,Di,L)
    alpha_i = Wärmeübergang(Nu,lamW,Di)
    kA = Wärmeduchgang(alpha_i,Di,Dm,Da,lamRohr1,lamRohr2)
    return kA
end

function Wärmeduchgang_Rohr_innen(kante)
    m = kante.y.m
    Di = kante.Param.Di
    Dm = kante.Param.Dm
    Da = kante.Param.Da
    A = kante.Param.A
    L = kante.Param.L
    mu = kante.Param.mu
    cv_H2O = kante.Param.cv_H2O
    lamW = kante.Param.lamW
    lamRohr1 = kante.Param.lamRohr1
    lamRohr2 = kante.Param.lamRohr2
    alpha_a = kante.Param.alpha_a
    Re = Reynolds(m,Di,A,mu)
    Pr = Prandtl(cv_H2O,mu,lamW)
    Nu = Nusselt(Re,Pr,Di,L)
    alpha_i = Wärmeübergang(Nu,lamW,Di)
    kA = Wärmeduchgang(alpha_i,alpha_a,Di,Dm,Da,lamRohr1,lamRohr2)
    return kA
end
#-----------------------

#-- Berechnung Wärmedurchgangskoeffizienten Wärmetauscher
function Wärmeduchgang_Wärmetauscher_Ringspalt(kante)
    m = kante.y.m
    Di = kante.Param.Di
    Dm = kante.Param.Dm
    Da = kante.Param.Da
    Dgl = Da-Dm
    A = kante.Param.Aa
    L = kante.Param.L
    mu = kante.Param.mu
    cv_H2O = kante.Param.cv_H2O
    lamW = kante.Param.lamW
    lamRohr = kante.Param.lamRohr
    Re = Reynolds(m,Dgl,A,mu)
    Pr = Prandtl(cv_H2O,mu,lamW)
    Nu = Nusselt(Re,Pr,Dgl,L)
    alpha_a = Wärmeübergang(Nu,lamW,Dgl)
    alpha_i = Wärmeübergang_Wärmetauscher_Rohrinnen(kante.R)
    kA = Wärmeduchgang(alpha_i,alpha_a,Di,Dm,lamRohr)
    kante.R.Z["kA"] = kA
    #@show kA
    return kA
end

function Wärmeübergang_Wärmetauscher_Rohrinnen(kante)
    m = kante.y.m
    D = kante.Param.Di
    A = kante.Param.Ai
    L = kante.Param.L
    mu = kante.Param.mu
    cv_H2O = kante.Param.cv_H2O
    lamW = kante.Param.lamW
    Re = Reynolds(m,D,A,mu)
    Pr = Prandtl(cv_H2O,mu,lamW)
    Nu = Nusselt(Re,Pr,D,L)
    alpha_i = Wärmeübergang(Nu,lamW,D)
    return alpha_i
end
#-----------------------

#-- Berechnung Wärmedurchgangskoeffizienten Wärmeübertragung
function Wärmeduchgang_Wärmeübertragung(kante)
    m = kante.y.m
    Di = kante.Param.Di
    Da = kante.Param.Da
    A = kante.Param.A
    L = kante.Param.L
    mu = kante.Param.mu
    cv_H2O = kante.Param.cv_H2O
    lamW = kante.Param.lamW
    lamRohr = kante.Param.lamRohr1
    Re = Reynolds(m,Di,A,mu)
    Pr = Prandtl(cv_H2O,mu,lamW)
    Nu = Nusselt(Re,Pr,Di,L)
    alpha_i = Wärmeübergang(Nu,lamW,Di)
    alpha_a = kante.Param.alpha_a
    kA = Wärmeduchgang(alpha_i,alpha_a,Di,Da,lamRohr)
    return kA
end
#-----------------------