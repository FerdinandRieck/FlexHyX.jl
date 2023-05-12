abstract type FlexHyX end

abstract type Knoten <: FlexHyX end
abstract type Kante <: FlexHyX end

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

function tuple2array(y_tuple)
    y_arr = Float64[]; P_scale = Float64[]; y_leg = String[]; 
    idx_ele = Dict()
    idx_ifluss = Array{Int}(undef, 0,2); 
    idx_mfluss = Array{Int}(undef, 0,2); 
    idx_efluss = Array{Int}(undef, 0,2); 
    k = 0;
    for i=1:length(y_tuple.knoten)
        for ff in fieldnames(typeof(y_tuple.knoten[i].y))
            if ff != :Param
                append!(y_arr,getfield(y_tuple.knoten[i].y,ff)); k +=1 
                if first(string(ff))=='P' P_scale = [P_scale; 1.0e-5]
                else P_scale = [P_scale; 0.0] end 
                leg_i = string(y_tuple.knoten[i].Z["Nr"],ff)
                y_leg = [y_leg; leg_i]
                idx_ele[leg_i] = [i k]  #-- Dictionary
            end
        end
    end
    for i=1:length(y_tuple.kanten)
        for ff in fieldnames(typeof(y_tuple.kanten[i].y))
            if ff != :Param
                append!(y_arr,getfield(y_tuple.kanten[i].y,ff)); k +=1
                if first(string(ff))=='i' idx_ifluss = [idx_ifluss;[i k]] end  #-- Strom der Kante i_k steht in y an Stelle k
                if first(string(ff))=='m' idx_mfluss = [idx_mfluss;[i k]] end  
                if first(string(ff))=='e' idx_efluss = [idx_efluss;[i k]] end 
                if first(string(ff))=='P' P_scale = [P_scale; 1.0e-5]
                else P_scale = [P_scale; 0.0] end 
                leg_i = string(y_tuple.kanten[i].Z["Nr"],ff)
                y_leg = [y_leg; leg_i]
                idx_ele[leg_i] = [i k]
            end
        end
    end
    return y_arr, idx_ifluss, idx_mfluss, idx_efluss, P_scale, y_leg, idx_ele
end

function array2tuple!(y_tuple,y_arr)
    idx = 0
    for i=1:length(y_tuple.knoten)
        for ff in fieldnames(typeof(y_tuple.knoten[i].y))
            if ff != :Param
                idx += 1
                setfield!(y_tuple.knoten[i].y, ff, y_arr[idx])
            end
        end
    end
    for i=1:length(y_tuple.kanten)
        for ff in fieldnames(typeof(y_tuple.kanten[i].y))
            if ff != :Param
                idx += 1
                setfield!(y_tuple.kanten[i].y, ff, y_arr[idx])
            end
        end
    end
    nothing
end