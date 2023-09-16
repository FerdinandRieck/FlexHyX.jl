function plotNetz(dir::String)
	#-- knoten_infos und kanten_infos erzeugen
    J_cfg = JSON.parsefile(dir*"/Netzwerk/flexhyx.cfg")
	pfad = get(J_cfg,"Pfad","."); pfad=dir*"/"*pfad*"/";
    netzfile = pfad*get(J_cfg,"Netzwerkfile",0)
    J = JSON.parsefile(netzfile)
	K = get(J,"Knoten",0)
	knoten_infos = []; kanten_infos = []
    for kk in K
		if haskey(kk,"#")==false && haskey(kk,"#Nr")==false
			push!(knoten_infos,kk);
        end
    end
	K = get(J,"Kanten",0)
	for kk in K
		if haskey(kk,"#")==false && haskey(kk,"#Nr")==false
			push!(kanten_infos,kk);
	   	end
	end
	#-- Von/Nach = aktualisieren, RefKante aktualisieren
	n_n = size(knoten_infos)[1];   n_e = size(kanten_infos)[1]
	nr2kn = zeros(Int,n_n); nr2ka = zeros(Int,n_e)
	for i = 1:n_n
		nr2kn[i] = knoten_infos[i]["Nr"];
	end
	for i = 1:n_e
		nr2ka[i] = kanten_infos[i]["Nr"];
	end
	for i = 1:n_e
		kanten_infos[i]["VonNach"][1] = findall(x->x==kanten_infos[i]["VonNach"][1],nr2kn)[1];
		kanten_infos[i]["VonNach"][2] = findall(x->x==kanten_infos[i]["VonNach"][2],nr2kn)[1];
	end
	#-------------------------------------

	P = plot(legend = false, ticks = false, showaxis = false, aspect_ratio=:equal, size=(2000,1000))
	plot_kanten!(knoten_infos,kanten_infos)
	plot_knoten!(knoten_infos)
	return P
end


function plot_knoten!(knoten)
    # Knoten darstellen
    n_knoten = length(knoten)
    xmin = Inf; xmax = -Inf; ymin = xmin; ymax = xmax
    for i = 1:n_knoten
        x = knoten[i]["XY"][1]; y = knoten[i]["XY"][2]
        xmin = min(xmin, x); ymin = min(ymin, y)
        xmax = max(xmax, x); ymax = max(ymax, y)
    end
    dx = (xmax - xmin + 10) / 200
    dy = (ymax - ymin + 10) / 200

	plot!(xlim = (xmin - 10 * dx, xmax + 10 * dx), ylims = (ymin - 10 * dy, ymax + 10 * dy))

    for i = 1:n_knoten
        kk = knoten[i]
        x = knoten[i]["XY"][1]; y = knoten[i]["XY"][2]
		s = Symbol(kk["Typ"],"_Knoten"); obj = getfield(FlexHyX, s)
		color = :blue
		if obj <: Strom_Knoten
			color = :red
		elseif obj <: Gas_Knoten
			color = :green
		elseif obj <: Temp_Knoten
			color = :magenta
		elseif obj <: Wasser_Knoten
			color = :blue
		end
		plot!([x],[y], seriestype = :scatter, markercolor = :white, markerstrokewidth = 2, markerstrokecolor = color)
		tt = string(kk["Nr"])*" "*kk["Typ"]
		fntsz = 9
		n = length(tt)
		str0 = "█"^(n÷2+1)
		annotate!(x+dx,y+dy, text(str0, :left, :white, fntsz+1))
		annotate!(x+dx,y+dy, text(tt, :left, color, fntsz))
	end	
end

function plot_kanten!(knoten, kanten)
    #-- Kanten plotten
    A = zeros(length(knoten),length(knoten))
    for i in eachindex(kanten)
		kk = kanten[i]
		iv = kk["VonNach"][1]; in = kk["VonNach"][2]
		x0 = knoten[iv]["XY"][1]; y0 = knoten[iv]["XY"][2]
		x1 = knoten[in]["XY"][1]; y1 = knoten[in]["XY"][2]
		xv = x1 - x0; yv = y1 - y0
		A[iv, in] += 1; A[in, iv] += 1
		iq = A[iv, in] - 1
		dquer = iq * 0.15 * ((-1)^iq)
		xq = -dquer * yv; yq = dquer * xv
		x13 = x0 + xv/3 + xq; x23 = x0 + 2*xv / 3 + xq
		y13 = y0 + yv/3 + yq; y23 = y0 + 2*yv/3 + yq
		x12 = x13+0.75*(x23-x13); y12 = y13+0.75*(y23-y13);
		color = :blue
		if (kk["Typ"] == "mE") || (kk["Typ"] == "mBZ")
			color = :green
		else
			s = Symbol(kk["Typ"],"_kante"); obj = getfield(FlexHyX, s)
			if obj <: Strom_Kante
				color = :red
			elseif obj <: Gas_Kante
				color = :green
			elseif obj <: Temp_Kante
				color = :red
			elseif obj <: Gas_Strom_Kante
				color = :orange
			elseif obj <: Wasser_Kante
				color = :blue
			end
		end
		GR.setarrowsize(0.5)
		P = plot!([x0, x13, x23, x1], [y0, y13, y23, y1],
			arrow=Plots.Arrow(:closed, :head, 2.5, 2.0),
			markersize = 20,
			color=color
		)
		tt = string(kk["Nr"])*" "*kk["Typ"]
		fntsz = 9
		n = length(tt)
		str0 = "█"^(n÷2+1)
		annotate!(x12, y12, text(str0, :white, fntsz+1))
		annotate!(x12, y12, text(tt, fntsz))
    end
end
