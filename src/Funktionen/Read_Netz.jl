function read_netz(dir,netzfile,zwerte,zt,znamen)
	typ_kn = ["U","U0"]
	typ_ka = ["iB","iV","iPV"]
	J = JSON.parsefile(netzfile)
	eventfile = []
	if haskey(J,"Events")
		eventfile = get(J,"Events",[]); eventfile = dir*"/Events/"*eventfile[2:end]*".jl"
		println("Eventfile:",eventfile)
		include(eventfile)
	end
    knoten = [];  kanten = []
	K = get(J,"Knoten",0)
    for kk in K
		if haskey(kk,"#")==false && haskey(kk,"#Nr")==false
     	    typ = kk["Typ"];
            if (typ=="U0") & (haskey(kk,"Spannung")==true) kk["U0"] = kk["Spannung"]; end
			if (typ=="GP0") & (haskey(kk,"T0")==true) kk["T0"] = kk["T0"] + 273.15; end
			if (typ=="T0") & (haskey(kk,"T0")==true) kk["T0"] = kk["T0"] + 273.15; end
			if (typ=="TM") & (haskey(kk,"T0")==true) kk["T0"] = kk["T0"] + 273.15; end
		#	if (typ=="GP0")|(typ=="GPSP")|(typ=="GPMH")... # besser direkt so
			if typ=="GPSP"
				if haskey(kk,"P0")==true kk["P0"] = kk["P0"]*1.0e5 end
				if haskey(kk,"T0")==true kk["T0"] = kk["T0"] + 273.15 end
			end
			if typ=="GPMH"
				if haskey(kk,"P0")==true kk["P0"] = kk["P0"]*1.0e5 end
				if haskey(kk,"T0")==true kk["T0"] = kk["T0"] + 273.15 end
				if haskey(kk,"Theta0")==true kk["Θ0"] = kk["Theta0"] end
			end
			push!(knoten,kk);
        end
    end
	K = get(J,"Kanten",0)
	for kk in K
		if haskey(kk,"#")==false && haskey(kk,"#Nr")==false
			typ = kk["Typ"];
			if (typ=="iPV")  
				if (haskey(kk,"Temp")==true) kk["T_PV"] = kk["Temp"] + 273.15; end
				if (haskey(kk,"Strahlung")==true) kk["G"] = kk["Strahlung"] end
			end
			#if (typ=="iV")&(haskey(kk,"R")==false) kk["R"]=0.0; end # ??? Wofür wird R benötigt???
			if (typ=="iE")&(haskey(kk,"Zellen")==true) kk["n_Z"] = kk["Zellen"] end
			if (typ=="mMH")&(haskey(kk,"Theta0")==true) kk["Θ0"] = kk["Theta0"] end
			for (k, v) in kk
				if v[1]=='@'
					fcn = v[2:end];
					#file = "Netzwerk/"*fcn*".jl";
					#println(file);   include(file)
					kk[k] = getfield(Main, Symbol(fcn))
				end
			end
			if (haskey(kk,"Schaltzeit")==true) kk["Schaltdauer"] = 60.0 end
			if (haskey(kk,"Zeitreihe")==true)
				spalte = first(findall(x->x==kk["Zeitreihe"],znamen))
				kk["zt"] = zt
				kk["zwerte"] = zwerte[:,spalte]
				kk["interpol"] = 0;	#-- Default
				kk["ym"] = copy(kk["zwerte"]); kk["ym"][1] = 0; #-- Integrale
				for i = 2:length(kk["zwerte"])
					trapez = 0.5*(kk["zwerte"][i] + kk["zwerte"][i-1]) * (kk["zwerte"][i] - kk["zwerte"][i-1]);
					kk["ym"][i] = kk["ym"][i-1] + trapez;
				end
				if (haskey(kk,"T")==true)
					kk["interpol"] = min(kk["T"],2);
					if kk["T"] > 2
					   ty, y, ym = reduktion(kk["zt"],zwerte[:,spalte],kk["T"])
					   kk["zt"] = ty; kk["zwerte"] = y; kk["ym"] = ym;
					end
				end
			end
			push!(kanten,kk);
	   end
	end
	
#-- Zweiter Durchlauf: Von/Nach = aktualisieren, RefKante aktualisieren
	n_n = size(knoten)[1];   n_e = size(kanten)[1]
	nr2kn = zeros(Int,n_n); nr2ka = zeros(Int,n_e)
	for i = 1:n_n
		nr2kn[i] = knoten[i]["Nr"];
	end
	for i = 1:n_e
		nr2ka[i] = kanten[i]["Nr"];
	end
	for i = 1:n_e
		kanten[i]["VonNach"][1] = findall(x->x==kanten[i]["VonNach"][1],nr2kn)[1];
		kanten[i]["VonNach"][2] = findall(x->x==kanten[i]["VonNach"][2],nr2kn)[1];
	end
	#-- RefKanten Infos in Kante einfügen und RefKanten Index überschreiben
	idx_ka = Int[]
	for i = 1:n_e
		if haskey(kanten[i],"RefKante")
			i_ka = findall(x->x==kanten[i]["RefKante"],nr2ka)[1];
			kanten[i]["RefKante"] = i_ka;
			RefKante = kanten[i_ka];
			kanten[i]["RefKante"] = RefKante	
			idx_ka = [idx_ka; i_ka]
		end
	end
	deleteat!(kanten,idx_ka)
#--
    return eventfile, knoten, kanten
end