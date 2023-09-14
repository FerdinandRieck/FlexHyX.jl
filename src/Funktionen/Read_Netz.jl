function readNetz(dir,netzfile,zwerte,zt,znamen)
	J = JSON.parsefile(netzfile)
	eventfile = []
	if haskey(J,"Events")
		eventfile = get(J,"Events",[]); eventfile = dir*"/Events/"*eventfile[2:end]*".jl"
		println("Eventfile:",eventfile)
		io = open(eventfile,"r"); close(io)
		include(eventfile)
	end
    knoten = [];  kanten = []
	K = get(J,"Knoten",0)
    for kk in K
		if haskey(kk,"#")==false && haskey(kk,"#Nr")==false
			typ = kk["Typ"]; 
			if haskey(kk,"T0")==true kk["T0"] = kk["T0"] + 273.15 end 
			if haskey(kk,"P0")==true kk["P0"] = kk["P0"]*1.0e5 end 
            if (typ=="U0") & (haskey(kk,"Spannung")==true) kk["U0"] = kk["Spannung"] end 
			if typ=="GPMH"
				if haskey(kk,"Theta0")==true kk["Θ0"] = kk["Theta0"] end
			end
			for (k, v) in kk
				if v[1]=='@'
					fcn = v[2:end];
					kk[k] = getfield(FlexHyX, Symbol(fcn))
				end
			end
			if (haskey(kk,"Schaltzeit")==true) && (haskey(kk,"Schaltdauer")==false)
				kk["Schaltdauer"] = 0.1 
			end
			push!(knoten,kk);
        end
    end
	K = get(J,"Kanten",0)
	for kk in K
		if haskey(kk,"#")==false && haskey(kk,"#Nr")==false
			typ = kk["Typ"]; 
			if (haskey(kk,"T0")==true) kk["T0"] = kk["T0"] + 273.15 end 
			if haskey(kk,"P0")==true kk["P0"] = kk["P0"]*1.0e5 end 
			if (typ=="iPV")  
				if (haskey(kk,"Temp")==true) kk["T_PV"] = kk["Temp"] + 273.15; end
				if (haskey(kk,"Strahlung")==true) kk["G"] = kk["Strahlung"] end
			end
			if (typ=="iE")&(haskey(kk,"Zellen")==true) kk["n_Z"] = kk["Zellen"] end
			if (typ=="mMH")&(haskey(kk,"Theta0")==true) kk["Θ0"] = kk["Theta0"] end
			for (k, v) in kk
				if v[1]=='@'
					fcn = v[2:end];
					kk[k] = getfield(FlexHyX, Symbol(fcn))
				end
			end
			if (haskey(kk,"Schaltzeit")==true) && (haskey(kk,"Schaltdauer")==false)
				kk["Schaltdauer"] = 0.1 
			end
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
	
#-- Von/Nach = aktualisieren, RefRohr aktualisieren
	n_n = size(knoten)[1];   n_e = size(kanten)[1]
	nr2kn = zeros(Int,n_n); nr2ka = zeros(Int,n_e)
	for i = 1:n_n
		nr2kn[i] = knoten[i]["Nr"];
	end
	for i = 1:n_e
		nr2ka[i] = kanten[i]["Nr"];
	end
	for i = 1:n_e
		if haskey(kanten[i],"RohrAustausch")
			i_ka = findall(x->x==kanten[i]["RohrAustausch"],nr2ka)[1];
			kanten[i]["RohrAustausch"] = i_ka 
		end
		if haskey(kanten[i],"Ringspalt")
			i_ka = findall(x->x==kanten[i]["Ringspalt"],nr2ka)[1];
			kanten[i]["Ringspalt"] = i_ka 
		end
		if haskey(kanten[i],"KnotenAussen")
			i_kn = findall(x->x==kanten[i]["KnotenAussen"],nr2kn)[1];
			kanten[i]["KnotenAussen"] = i_kn 
		end
		kanten[i]["VonNach"][1] = findall(x->x==kanten[i]["VonNach"][1],nr2kn)[1];
		kanten[i]["VonNach"][2] = findall(x->x==kanten[i]["VonNach"][2],nr2kn)[1];
	end
	#-- RefKanten Infos in Kante einfügen und RefKante löschen
	idx_ka = Int[]
	for i = 1:n_e
		if haskey(kanten[i],"RefKante")
			i_ka = findall(x->x==kanten[i]["RefKante"],nr2ka)[1];
			RefKante = kanten[i_ka];
			kanten[i]["RefKante"] = RefKante	
			push!(idx_ka,i_ka)
		end
	end
	deleteat!(kanten,idx_ka)
#--
    return knoten, kanten, eventfile
end