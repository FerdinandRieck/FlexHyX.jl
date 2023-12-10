function netzwerk2y(knoten,kanten)
    y_arr = Float64[]
    idx_ele = Dict()
    idx_iflussL = Array{Int}(undef, 0,2); idx_iflussR = Array{Int}(undef, 0,2); 
    idx_mflussR = Array{Int}(undef, 0,2); idx_mflussL = Array{Int}(undef, 0,2); 
    idx_eflussR = Array{Int}(undef, 0,2); idx_eflussL = Array{Int}(undef, 0,2);
    k = 0;
    for i in eachindex(knoten)
        for ff in fieldnames(typeof(knoten[i].y))
            if ff != :Param
                n = length(getfield(knoten[i].y,ff))
                k = k + n
                append!(y_arr,getfield(knoten[i].y,ff));              
                leg_i = string(ff,knoten[i].Z["Nr"])
                idx_ele[leg_i] = [i k]  #-- Dictionary
            end
        end
    end
    for i in eachindex(kanten)
        for ff in fieldnames(typeof(kanten[i].y))
            if ff != :Param
                n = length(getfield(kanten[i].y,ff))
                k = k + n
                append!(y_arr,getfield(kanten[i].y,ff)); 
                if first(string(ff))=='i'
                    if string(ff)=="iR"
                        idx_iflussR = vcat(idx_iflussR, [i k])
                    elseif string(ff)=="iL"
                        idx_iflussL = vcat(idx_iflussL, [i k])
                    else
                        idx_iflussR = vcat(idx_iflussR, [i k])
                        idx_iflussL = vcat(idx_iflussL, [i k])
                    end
                end
                if first(string(ff))=='m' 
                    if string(ff)=="mR"
                        idx_mflussR = vcat(idx_mflussR, [i k])
                    elseif string(ff)=="mL"
                        idx_mflussL = vcat(idx_mflussL, [i k])
                    else
                        idx_mflussR = vcat(idx_mflussR, [i k])
                        idx_mflussL = vcat(idx_mflussL, [i k])
                    end
                end
                if first(string(ff))=='e'
                    if string(ff)=="eR"
                        idx_eflussR = vcat(idx_eflussR, [i k])
                    elseif string(ff)=="eL"
                        idx_eflussL = vcat(idx_eflussL, [i k])
                    else
                        idx_eflussR = vcat(idx_eflussR, [i k])
                        idx_eflussL = vcat(idx_eflussL, [i k])
                    end
                end
                leg_i = string(kanten[i].Z["Nr"],ff)
                idx_ele[leg_i] = [i k]  #-- Dictionary
            end
        end
    end
    return y_arr, idx_iflussL, idx_iflussR, idx_mflussL, idx_mflussR, idx_eflussL, idx_eflussR, idx_ele
end

function y2netzwerk!(knoten,kanten,y_arr)
    idx = 1
    for i in eachindex(knoten)
        for ff in fieldnames(typeof(knoten[i].y))
            if ff != :Param
                setfield!(knoten[i].y, ff, y_arr[idx])
                idx += 1
            end
        end
    end
    for i in eachindex(kanten)
        for ff in fieldnames(typeof(kanten[i].y))
            if ff != :Param
                n = length(getfield(kanten[i].y,ff))
                if n == 1
                    setfield!(kanten[i].y, ff, y_arr[idx])
                else
                    setfield!(kanten[i].y, ff, y_arr[idx:idx+n-1])
                end
                idx = idx + n
            end
        end
    end
    nothing
end

#-- Wird für Leistung_anhängen benötigt
function idx2netzwerk!(knoten,kanten)
    idx = 1
    for i in eachindex(knoten)
        for ff in fieldnames(typeof(knoten[i].y))
            if ff != :Param
                setfield!(knoten[i].y, ff, idx)
                idx += 1
            end
        end
    end
    for i in eachindex(kanten)
        for ff in fieldnames(typeof(kanten[i].y))
            if ff != :Param
                n = length(getfield(kanten[i].y,ff))
                if n == 1
                    setfield!(kanten[i].y, ff, idx)
                else
                    setfield!(kanten[i].y, ff, (idx:idx+n-1))
                end
                idx = idx + n
            end
        end
    end
    nothing
end