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
		x = x+dx*rand(2); y = y+dy*rand(2)
		P=plot!(x,y,line=:arrow,label=kanten[i]["Typ"])
	end
#	display(P)
	return P
end