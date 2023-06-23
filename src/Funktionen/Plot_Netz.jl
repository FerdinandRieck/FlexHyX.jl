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

	P = plot(legend = false, ticks = false, showaxis = false, aspect_ratio=:equal)
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
		if obj <: Strom_Knoten
			color = :red
		elseif obj <: Gas_Knoten
			color = :blue
		elseif obj <: Temp_Knoten
			color = :magenta
		end
		plot!([x],[y], seriestype = :scatter, markercolor = :white, markerstrokewidth = 2, markerstrokecolor = color)
		tt = string(kk["Nr"])*" "*kk["Typ"]
		fntsz = 6
		n = length(tt)
		str0 = "█"^(n÷2+1)
		annotate!(x+dx,y+dy, text(str0, :left, :white, fntsz+1))
		annotate!(x+dx,y+dy, text(tt, :left, color, fntsz))
	end	
end

function plot_kanten!(knoten, kanten)
    #-- Kanten plotten
    A = zeros(length(knoten),length(knoten))
    for i = 1:length(kanten)
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
		if (kk["Typ"] == "mE") || (kk["Typ"] == "mBZ")
			color = :green
		else
			s = Symbol(kk["Typ"],"_kante"); obj = getfield(FlexHyX, s)
			if obj <: Strom_Kante
				color = :blue
			elseif obj <: Gas_Kante
				color = :green
			elseif obj <: Temp_Kante
				color = :red
			elseif obj <: Gas_Strom_Kante
				color = :orange
			end
		end
		P = plot!([x0, x13, x23, x1], [y0, y13, y23, y1],
			arrow=Plots.Arrow(:closed, :head, 2.5, 2.0),
			color=color
		)
		tt = string(kk["Nr"])*" "*kk["Typ"]
		fntsz = 6
		n = length(tt)
		str0 = "█"^(n÷2+1)
		annotate!(x12, y12, text(str0, :white, fntsz+1))
		annotate!(x12, y12, text(tt, fntsz))
    end
end

#=
function plotNetz(knoten,kanten)
    n_n = size(knoten)[1]; n_e = size(kanten)[1]; iart_max = 0;
    for i=1:n_n
	  iart = -knoten[i]["iart"]; iart_max=max(iart_max,iart);
    end
    n_art = zeros(Int,iart_max); t_art = Array{String}(undef,iart_max)
    x = zeros(iart_max,n_n); y=zeros(iart_max,n_n)
    xmin=Inf; xmax=-Inf; ymin=Inf; ymax=-Inf;
    for i=1:n_n
        xy = knoten[i]["XY"]; iart = -knoten[i]["iart"]
	    n_art[iart] = n_art[iart]+1
        x[iart,n_art[iart]] = xy[1]; y[iart,n_art[iart]] = xy[2]
	    t_art[iart] = knoten[i]["Typ"]
	    xmin=min(xmin,xy[1]); xmax=max(xmax,xy[1])
	    ymin=min(ymin,xy[2]); ymax=max(ymax,xy[2])
    end
	P = plot()
	for i=1:iart_max
	    if n_art[i]>0
	        xx = x[i,1:n_art[i]];  yy = y[i,1:n_art[i]];
            P=plot!(xx,yy,seriestype=:scatter,title="FlexHyX",label=t_art[i])
        end
	end
 	x = zeros(2); y = zeros(2)
	dx = 0.01*(xmax-xmin); dy = 0.01*(ymax-ymin)
	for i=1:n_e
	    i0 = kanten[i]["VonNach"][1]; i1 = kanten[i]["VonNach"][2]
		x[1] = knoten[i0]["XY"][1]; x[2] = knoten[i1]["XY"][1]
		y[1] = knoten[i0]["XY"][2]; y[2] = knoten[i1]["XY"][2]
		#x = x+dx*rand(2); y = y+dy*rand(2)
		GR.setarrowsize(1)
		P=plot!(x,y,
			#arrow=arrow(:closed, 2.0),
			label=kanten[i]["Typ"],
			arrow=Plots.Arrow(:closed, :both, 2.5, 2.0),
			#bar_edges = true
			)
	end
	display(P)
end
=#

#=
function plotNetz(dir::String)
	#-- knoten_infos und kanten_infos erzeugen
    J_cfg = JSON.parsefile(dir*"/Netzwerk/flexhyx.cfg")
    now = Dates.now(); jetzt = [Dates.year(now) Dates.month(now) Dates.day(now) Dates.hour(now) Dates.minute(now) 0]
    startzeit = get(J_cfg,"Startzeit",jetzt)
    startzeit = String(Symbol(startzeit'))
    startzeit = startzeit[5:end-1]
    startzeit = DateTime(startzeit,"yyyy mm dd HH MM SS")
    pfad = get(J_cfg,"Pfad","."); pfad=dir*"/"*pfad*"/";
    netzfile = pfad*get(J_cfg,"Netzwerkfile",0)
    zeitfile = get(J_cfg,"Zeitreihenfile",nothing) 

    znamen = []; zwerte = []; zt = [];
    if zeitfile != nothing
        zstart, zt, zwerte, znamen, zeinheit, ztitel = readZeitreihe(pfad*zeitfile)
        dt = Second(zstart-startzeit); dt = dt.value
        zt = zt .+ dt
    end

    knoten_infos, kanten_infos, eventfile = readNetz(dir, netzfile, zwerte, zt, znamen)
	#------------------------------------------

	g = DiGraph(length(knoten_infos),0); X = []; Y = []; name = String[]
	for i = 1:length(knoten_infos)
		kk = knoten_infos[i]
		x = kk["XY"][1]; push!(X,x)
		y = kk["XY"][2]; push!(Y,y)
		typ = rpad(kk["Typ"],4)
		push!(name,typ)
	end
	#X = 1.7 * X; Y = 1.7 * Y
	for i = 1:length(kanten_infos)
		kk = kanten_infos[i]
		von = kk["VonNach"][1]
		nach = kk["VonNach"][2]
		add_edge!(g,von,nach)
		if (kk["Typ"] == "iE") || (kk["Typ"] == "iBZ")
			RefKante = kk["RefKante"]
			von_Ref = RefKante["VonNach"][1]; nach_Ref = RefKante["VonNach"][2]
			add_edge!(g,von_Ref,nach_Ref)
		end
	end
	graphplot(g,
	#size=(1000,1000),
    #-- Knoten
    names = name,
    fontsize = 7,
    nodeshape = :rect,
    markersize = 0.2,
    x = X, y = Y,
    #-- Kanten
    #edgelabel = edgelabels,
    curves = false,
	curvature_scalar = 0.05,
    )
	#display(P)
end
=#