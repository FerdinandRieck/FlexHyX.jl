#function read_zeitreihe(file::AbstractString)
function readZeitreihe(file)
    # Daten einer/mehrerer Zeitreihen einlesen
    T = read(file, String)
    J = JSON.parse(T)
    format = J["Format"]
    dt = 0.0
    if format == "Const_Delta"
        dt = J["dt"]
    end
    titel = J["Titel"]
    namen = J["Namen"]
    einheit = get(J, "Einheiten", [])
    startzeit = get(J, "Startzeit", nothing)
    if startzeit != nothing
        startzeit = String(Symbol(startzeit'))
        startzeit = startzeit[5:end-1]
        startzeit = DateTime(startzeit,"yyyy mm dd HH MM SS")
    end
    n = length(J["Daten"])
    t = zeros(n)
    werte = zeros(n, length(namen)) 
    for i = 1:n
        if format == "JMTSMS"
            date = J["Daten"][i]["JMTSMS"]
            date = String(Symbol(date'))
            date = date[5:end-1]
            if startzeit == nothing
                startzeit = DateTime(date,"yyyy mm dd HH MM SS")
            end
            dt = Second(DateTime(date,"yyyy mm dd HH MM SS") - startzeit)
            t[i] = dt.value
        end
        if format == "ISO"
            date = J["Daten"][i]["ISO"]
            if startzeit == nothing
                startzeit = Dates.DateTime(date)
            end
            t[i] = Dates.value(Dates.Time(date, "yyyy-mm-ddTHH:MM:SS.sssZ")) - Dates.value(Dates.Time(startzeit))
        end
        if format == "Delta"
            t[i] = J["Daten"][i]["Delta"]
        end
        if format == "Const_Delta"
            t[i] = (i - 1) * dt
        end
        JW = J["Daten"][i]["Wert"]
        # ??? allownan=true gibt es auch als option für JSON.parse ???
        if all(isa.(JW,Number))
            werte[i, :] = JW'
        else #-- NaN als Text vorhanden
            ff = JW .== "NaN"
            werte[i,!ff] = parse.(Float64, JW[.!ff])
            werte[i,ff] = NaN
        end
    end 
    #--
    if any(isnan.(werte))
        for i = 1:length(namen)
            werte[:,i] = fillmissing(werte[:,i], "linear", SortedDict(t => 1:length(t)))
        end
    end
    startdatum = startzeit
    #--
    return startdatum, t, werte, namen, einheit, titel
end

function getwert(t,zt,zwerte,i_art,y_m)
    #-- i_art = 0: linear, 1:pchip, 2:sigmoid, -1:werte y_m aus
    #--
    n = length(zt)
    dt = zt[2] - zt[1]
    lasti = round(Int,(t - zt[1]) / dt)
    lasti = max(min(lasti, n), 1)
    if t <= zt[1]
        w = zwerte[1]
        if i_art == -1
            w = y_m[1]
        end
    elseif t >= zt[n]
        w = zwerte[n]
        if i_art == -1
            w = y_m[n]
        end
    else
        if zt[lasti] > t
            while (lasti > 1) && (zt[lasti] > t)
                lasti = lasti - 1
            end
        else
            while (lasti < n - 1) && (zt[lasti+1] <= t)
                lasti = lasti + 1
            end
        end
        t0 = zt[lasti]
        t1 = zt[lasti+1]
        dt = t1 - t0
        y0 = zwerte[lasti]
        y1 = zwerte[lasti+1]
        if i_art == 0
            w = y0 + (t-t0)/dt * (y1-y0)
        elseif i_art == -1
            w = y_m[lasti] + (t-t0)/dt * (y_m[lasti+1]-y_m[lasti])
        elseif i_art == 1
            l0 = max(1, lasti-1)
            l1 = min(lasti+2, n)
            pp = pchip(zt[l0:l1], zwerte[l0:l1])
            w = pp(t)
        elseif i_art == 2
            #x0 = 0.5
            #if exists("ym") && !isempty(ym)
            x0 = 1 - ((y_m[lasti+1]-y_m[lasti])/dt-y0)/(y1-y0)
            x0 = min(max(x0, 0.1), 0.9)
            #end
            w = y0 + (y1-y0)*sig((t-t0)/dt, x0)
        end
    end
    #--
    return w
end

function sig(x,x0)
    # ep = 1.0e-5; b = log(1/ep-1);
    b = 11.512915464920228
    a = b/min(x0,1-x0); y = 1 ./ (1 + exp(-a*(x-x0)))
    # a = atanh(1-ep); s = a/min(x0,1-x0); y = (tanh(s*(x-x0))+1)/2;
    # dauer = min(x0,1-x0); u = (x-x0+dauer)/(2*dauer); u=min(max(u,0),1); y = 3*u.^2 - 2*u.^3; 
    return y
end

# function reduktion(t::Vector{Typ}, x::Vector{Typ}, T::Typ) where {Typ <: Real}
function reduktion(t, x, T) 
    #-- Datenreduktion
    t0 = t[1]; t1 = t0 + T; n = length(t)
    k = 1; y_min = typemax(T); y_max = -typemax(T)
    ind = similar(1:n); m_n = similar(1:n)
    for i=1:n
        if x[i] < y_min
            ind[k] = i
            y_min = x[i]
            m_n[k] = -1
        end
        if x[i] >= y_max
            ind[k+1] = i
            y_max = x[i]
            m_n[k+1] = 1
        end
        # !!! Unterschiede bei Matlab und Julia mit >=
        if t[i] >= t1 || i == n #-- Zeitraum fertig
            k += 2
            t1 += T
            y_min = typemax(T)
            y_max = -typemax(T)
        end
    end
    ind = ind[1:k-1]
    for i=1:minimum(ind)-1
        ind = [ind; i]
        m_n = [m_n; 0]
        #k += 1
    end
    for i=maximum(ind)+1:n
        ind = [ind; i]
        m_n = [m_n; 0]
        #k += 1
    end
   # max_k = k-1
    ty = t[ind]; y = x[ind]
    isort = sortperm(ty)
    ty = ty[isort]; y = y[isort]; m_n = m_n[isort]; ind = ind[isort]
    #-- bei Übergang min-min oder max-max Werte einfügen
    for k = 1:length(ty)-1
        if m_n[k]*m_n[k+1] > 0
            y_min = typemax(T)
            y_max = -typemax(T)
            i_min = 0
            i_max = 0
            for i = ind[k]+1:ind[k+1]-1
                if x[i] < y_min
                    y_min = x[i]
                    i_min = i
                end
                if x[i] > y_max
                    y_max = x[i]
                    i_max = i
                end
            end
            if m_n[k] > 0 && i_min > 0
               # max_k += 1
               # ind[max_k] = i_min
               ind = [ind; i_min]
            end
            if m_n[k] < 0 && i_max > 0
              #  max_k += 1
              #  ind[max_k] = i_max # !!! prüfen ob die Position von max_k relevant ist
                ind = [ind; i_max]
            end
        end
    end
    ty = t[ind]; y = x[ind]
    isort = sortperm(ty)
    ty = ty[isort]; y = y[isort]
    #-- Integrale zwischen den Extrema berechnen
    k = 1; sum = 0.0; ym = zeros(n); ym[1] = 0.0
    for i = 2:n
        trapez = 0.5*(x[i]+x[i-1])*(t[i]-t[i-1])
        sum += trapez
        if ty[k+1] == t[i] #-- neuer Zeitraum
            ym[k+1] = ym[k] + sum
            k += 1
            sum = 0.0
        end
    end
    ym = ym[1:k]
    return ty, y, ym
end
