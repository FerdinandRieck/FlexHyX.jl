#-- PtGtX
function plotSol(y,t)
  n_par = 1
  p1 = plot(t/3600,y[50,:], linewidth = 2, xlabel = "Zeit /h", ylabel = "W", title = "Leistung", label = "2iV-P",legend=:topright, legend_columns=-1)
  p1 = plot!(t/3600,y[51,:], linewidth = 2, label = "3iPV-P")
  p2 = plot(t/3600,y[26,:], linewidth = 2, xlabel = "Zeit /h", ylabel = "A", title = "Ströme", label = "1iB",legend=:bottom, legend_columns=-1)
  p2 = plot!(t/3600,y[34,:], linewidth = 2, label = "6iE")
  p2 = plot!(t/3600,y[39,:], linewidth = 2, label = "9iBZ")
  p4 = plot(t/3600,y[28,:]/(18.3*3600*n_par), linewidth = 2, xlabel = "Zeit /h", ylabel = "%", title = "Ladezustand", label = "soc",legend=:bottom)
  p5 = plot(t/3600,y[22,:], linewidth = 2, xlabel = "Zeit /h", ylabel = "%", title = "Beladung MHS", label = "Θ",legend=:top)

  p = plot(p1,p2,p4,p5, layout = (4, 1), size=(1000,1000),  legendfontsize= 12)
  return p
end

#-- Fernwärme
function plotSol(y,t,idx_ele,Typ)
  if Typ == "Zweileiternetzwerk"
    T1 = idx_ele["T1"][2] #-- Vorlauf Fern
    T_5 = idx_ele["5T"][2]  #-- Rücklauf Fern
    T9 = idx_ele["T9"][2] #-- Raumtemp. 1
    T6 = idx_ele["T6"][2]   #-- Vorlauf Haus 1
    T_10 = idx_ele["10T"][2] #-- Rücklauf Haus 1
    T_a = idx_ele["T10"][2]


    p1 = plot(t/(3600*24),y[T1,:].-273.15, linewidth = 1.5, label = "Vorlauf Fernwärme")
    plot!(t/(3600*24),y[T6,:].-273.15, linewidth = 1.5, label = "Vorlauf Haus 1")
    plot!(t/(3600*24),y[T9,:].-273.15, linewidth = 1.5, label = "Raumtemp. Haus 1")
    plot!(t/(3600*24),y[T_5,:].-273.15, linewidth = 1.5, label = "Rücklauf Fernwärme")
    plot!(t/(3600*24),y[T_10,:].-273.15, linewidth = 1.5, label = "Rücklauf Haus 1")
    plot!(t/(3600*24),y[T_a,:].-273.15, linewidth = 1.5, label = "Außentemp.")

    p = plot(p1, 
    layout = (1,1),
    left_margin=1mm, 
    size=(1000*0.7,800*0.7),  
    xlabel = "Zeit /Tage", ylabel = "°C", title = "Temperaturen", 
    xlims=(0,7), 
    ylims=(-20,70), 
    legend=:outerbottom,
    legend_columns=3,
    legendfontsize= 10,
    tickfont= 9)
  end
  if Typ == "Einleiternetzwerk"
    T1 = idx_ele["T1"][2] #-- Vorlauf Fern
    T_4 = idx_ele["4T"][2]  #-- Rücklauf Fern
    T8 = idx_ele["T8"][2] #-- Raumtemp. 1
    T5 = idx_ele["T5"][2]   #-- Vorlauf Haus 1
    T11 = idx_ele["T11"][2] #-- Rücklauf Haus 1
    T_a = idx_ele["T9"][2]


    p1 = plot(t/(3600*24),y[T1,:].-273.15, linewidth = 1.5, label = "Vorlauf Fernwärme")
    plot!(t/(3600*24),y[T5,:].-273.15, linewidth = 1.5, label = "Vorlauf Haus 1")
    plot!(t/(3600*24),y[T8,:].-273.15, linewidth = 1.5, label = "Raumtemp. Haus 1")
    plot!(t/(3600*24),y[T_4,:].-273.15, linewidth = 1.5, label = "Rücklauf Fernwärme")
    plot!(t/(3600*24),y[T11,:].-273.15, linewidth = 1.5, label = "Rücklauf Haus 1")
    plot!(t/(3600*24),y[T_a,:].-273.15, linewidth = 1.5, label = "Außentemp.")

    p = plot(p1, 
    layout = (1,1),
    left_margin=1mm, 
    size=(1000*0.7,800*0.7),  
    xlabel = "Zeit /Tage", ylabel = "°C", title = "Temperaturen", 
    xlims=(0,7), 
    ylims=(-20,70), 
    legend=:outerbottom,
    legend_columns=3,
    legendfontsize= 10,
    tickfont= 9)
  end
  return p
end