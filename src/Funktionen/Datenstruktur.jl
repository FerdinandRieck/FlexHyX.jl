abstract type flexhyx end

abstract type Knoten <: flexhyx end
abstract type Kante <: flexhyx end

abstract type Strom_Knoten <: Knoten end
abstract type Gas_Knoten <: Knoten end
abstract type Temp_Knoten <: Knoten end

abstract type Strom_Kante <: Kante end
abstract type Gas_Kante <: Kante end
abstract type Temp_Kante <: Kante end
abstract type Gas_Strom_Kante <: Kante end

Base.@kwdef mutable struct Netzwerk 
    knoten
    kanten
end

function netzwerk2array(y_netz)
    y_arr = Float64[]; P_scale = Float64[]; 
    leg_knoten = String[]; leg_kanten = String[]; 
    idx_ele = Dict()
    idx_iflussL = Array{Int}(undef, 0,2);
    idx_iflussR = Array{Int}(undef, 0,2); 
    idx_mfluss = Array{Int}(undef, 0,2); 
    idx_efluss = Array{Int}(undef, 0,2); 
    k = 0;
    for i=1:length(y_netz.knoten)
        for ff in fieldnames(typeof(y_netz.knoten[i].y))
            if ff != :Param
                append!(y_arr,getfield(y_netz.knoten[i].y,ff)); k +=1 
                if first(string(ff))=='P' push!(P_scale, 1.0e-5)
                else push!(P_scale, 0.0) end             
                leg_i = string(y_netz.knoten[i].Z["Nr"],ff)
                push!(leg_knoten, leg_i*"-"*y_netz.knoten[i].Z["Typ"])
                idx_ele[leg_i] = [i k]  #-- Dictionary
            end
        end
    end
    for i=1:length(y_netz.kanten)
        for ff in fieldnames(typeof(y_netz.kanten[i].y))
            if ff != :Param
                append!(y_arr,getfield(y_netz.kanten[i].y,ff)); k +=1
                if first(string(ff))=='i'
                    idx_iflussL = vcat(idx_iflussL, [i k]) #-- Strom der Kante i_k steht in y an Stelle k
                    idx_iflussR = vcat(idx_iflussR, [i k])
                    if string(ff)=="i_out"
                        idx_iflussL = idx_iflussL[1:end-1,:]
                        idx_iflussR = idx_iflussR[1:end .!= end-1,:] #-- Lösche vorletzte Zeile
                    end
                end
                if first(string(ff))=='m' idx_mfluss = vcat(idx_mfluss, [i k]) end  
                if first(string(ff))=='e' idx_efluss = vcat(idx_efluss, [i k]) end 
                if first(string(ff))=='P' push!(P_scale, 1.0e-5)
                else push!(P_scale, 0.0) end 
                leg_i = string(y_netz.kanten[i].Z["Nr"],ff)
                push!(leg_kanten, leg_i*"-"*y_netz.kanten[i].Z["Typ"])
                idx_ele[leg_i] = [i k]  #-- Dictionary
            end
        end
    end
    return y_arr, idx_iflussL, idx_iflussR, idx_mfluss, idx_efluss, P_scale, leg_knoten, leg_kanten, idx_ele
end

function array2netzwerk!(y_netz,y_arr)
    idx = 0
    for i=1:length(y_netz.knoten)
        for ff in fieldnames(typeof(y_netz.knoten[i].y))
            if ff != :Param
                idx += 1
                setfield!(y_netz.knoten[i].y, ff, y_arr[idx])
            end
        end
    end
    for i=1:length(y_netz.kanten)
        for ff in fieldnames(typeof(y_netz.kanten[i].y))
            if ff != :Param
                idx += 1
                setfield!(y_netz.kanten[i].y, ff, y_arr[idx])
            end
        end
    end
    nothing
end

function idx2netzwerk!(y_netz)
    idx = 0
    for i=1:length(y_netz.knoten)
        for ff in fieldnames(typeof(y_netz.knoten[i].y))
            if ff != :Param
                idx += 1
                setfield!(y_netz.knoten[i].y, ff, idx)
            end
        end
    end
    for i=1:length(y_netz.kanten)
        for ff in fieldnames(typeof(y_netz.kanten[i].y))
            if ff != :Param
                idx += 1
                setfield!(y_netz.kanten[i].y, ff, idx)
            end
        end
    end
    nothing
end