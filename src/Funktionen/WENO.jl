#-- WENO alt
function recover!(y,yL,yR)
    if length(y) > 1
        yL[1] = 2*y[1]-y[2];  yL[2:end] = y
        yR[end]= 2*y[end]-y[end-1]; yR[1:end-1] = y
    end
    return nothing
end

function recover_weno!(y,yL,yR)
    #-- Zellenmittelwerte auf Zellgrenzen interpolieren
    #-- L = upwind, R = downwind, WENO 3. Ordnung
    n = length(y);
    ep = 1.0e-6; p = 0.6;
    if n > 1
        yL[1] = 11/6*y[1]-7/6*y[2]+y[3]/3; yR[1] = yL[1]; #-- Randwerte
        yL[2] = y[1]/3+5/6*y[2]-y[3]/6;
        yL[n+1] = 11/6*y[n]-7/6*y[n-1]+y[n-2]/3; yR[n+1] = yL[n+1];
        yR[n] = y[n]/3+5/6*y[n-1]-y[n-2]/6; 
        for i=2:n-1
            uL = y[i]-y[i-1]; uC = y[i+1]-2*y[i]+y[i-1]; uR = y[i+1]-y[i]; uCC = y[i+1]-y[i-1];
            ISL = uL^2; ISC = 13/3*uC^2 +0.25*uCC^2; ISR = uR^2; aL = 0.25*(1/(ep+ISL))^p; aC = 0.5*(1/(ep+ISC))^p;
            aR = 0.25*(1/(ep+ISR))^p;
            suma = max(aL+aC+aR,eps(1.0)); wL = aL/suma; wC = aC/suma; wR = aR/suma;
            yR[i] = (0.5*wL+5/12*wC)*y[i-1] + (0.5*wL+2/3*wC+1.5*wR)*y[i] + (-wC/12-0.5*wR)*y[i+1];
            yL[i+1] = (-0.5*wL-wC/12)*y[i-1] + (1.5*wL+2/3*wC+0.5*wR)*y[i] + (5/12*wC+0.5*wR)*y[i+1];
        end 
    end
    return nothing
end
#---

#-- WENO mit LaxFriedrich ansatz
function recover(y)
    if length(y) > 1
        yL = [2*y[1]-y[2]; y]; yR = [y; 2*y[end]-y[end-1]] 
    else
        yL = [y[1]; y[1]]; yR = [y[end]; y[end]]
    end
    return yL, yR
end

function recover_weno(y)
    #-- Zellenmittelwerte auf Zellgrenzen interpolieren
    #-- L = upwind, R = downwind, WENO 3. Ordnung
    n = length(y);
    yL = Array{Number}(undef, n+1); yR = Array{Number}(undef, n+1) 
    #yL[1] = 11/6*y[1]-7/6*y[2]+y[3]/3; yR[1] = yL[1]; #-- Randwerte
    #yL[2] = y[1]/3+5/6*y[2]-y[3]/6;
    #yL[n+1] = 11/6*y[n]-7/6*y[n-1]+y[n-2]/3; yR[n+1] = yL[n+1];
    #yR[n] = y[n]/3+5/6*y[n-1]-y[n-2]/6; 
    yL[1] = 2*y[1]-y[2]; yR[1] = yL[1]; yL[2] = y[1]
    yR[n+1] = 2*y[n]-y[n-1]; yL[n+1] = yR[n+1]; yR[n] = y[n] 
    for i=2:n-1
        yR[i], yL[i+1] = weno3(y[i-1:i+1]); 
    end
    return yL, yR 
end

function weno3(y) #-- y = [y1,y2,y3]
    ep = 1.0e-6; p = 0.6;
    uL = y[2]-y[1]; uC = y[3]-2*y[2]+y[1]; uR = y[3]-y[2]; uCC = y[3]-y[1];
    ISL = uL^2; ISC = 13/3*uC^2 +0.25*uCC^2; ISR = uR^2; 
    aL = 0.25*(1/(ep+ISL))^p; aC = 0.5*(1/(ep+ISC))^p;  aR = 0.25*(1/(ep+ISR))^p;
    suma = max(aL+aC+aR,eps(1.0)); wL = aL/suma; wC = aC/suma; wR = aR/suma;
    y12 = (0.5*wL+5/12*wC)*y[1] + (0.5*wL+2/3*wC+1.5*wR)*y[2] + (-wC/12-0.5*wR)*y[3];
    y23 = (-0.5*wL-wC/12)*y[1] + (1.5*wL+2/3*wC+0.5*wR)*y[2] + (5/12*wC+0.5*wR)*y[3];
    return y12, y23
end

