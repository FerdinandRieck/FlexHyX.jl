function sum_i(IN,OUT)
    sum_IN = 0.0; sum_OUT = 0.0
    for i in eachindex(IN)
        if isempty(IN) == false
            if hasfield(typeof(IN[i].y), :iR)
                sum_IN = sum_IN + IN[i].y.iR
            elseif hasfield(typeof(IN[i].y), :i)
                sum_IN = sum_IN + IN[i].y.i
            end
        end
    end
    for i in eachindex(OUT)
        if isempty(OUT) == false
            if hasfield(typeof(OUT[i].y), :iL)
                sum_OUT = sum_OUT + OUT[i].y.iL
            elseif hasfield(typeof(OUT[i].y), :i)
                sum_OUT = sum_OUT + OUT[i].y.i
            end
        end
    end
    sum_i = sum_IN - sum_OUT
end

function sum_m(IN,OUT)
    sum_IN = 0.0; sum_OUT = 0.0
    for i in eachindex(IN)
        if isempty(IN) == false
            if hasfield(typeof(IN[i].y), :mR)
                sum_IN = sum_IN + IN[i].y.mR
            elseif hasfield(typeof(IN[i].y), :m)
                sum_IN = sum_IN + IN[i].y.m
            end
        end
    end
    for i in eachindex(OUT)
        if isempty(OUT) == false
            if hasfield(typeof(OUT[i].y), :mL)
                sum_OUT = sum_OUT + OUT[i].y.mL
            elseif hasfield(typeof(OUT[i].y), :m)
                sum_OUT = sum_OUT + OUT[i].y.m
            end
        end
    end
    sum_m = sum_IN - sum_OUT
end

function sum_e(IN,OUT)
    sum_IN = 0.0; sum_OUT = 0.0
    for i in eachindex(IN)
        if isempty(IN) == false
            if hasfield(typeof(IN[i].y), :eR)
                sum_IN = sum_IN + IN[i].y.eR
            elseif hasfield(typeof(IN[i].y), :e)
                sum_IN = sum_IN + IN[i].y.e
            end
            if haskey(IN[i].Z,"eR")
                sum_IN = sum_IN + IN[i].Z["eR"]
            end
            if haskey(IN[i].Z,"e")
                sum_IN = sum_IN + IN[i].Z["e"]
            end
        end
    end
    for i in eachindex(OUT)
        if isempty(OUT) == false
            if hasfield(typeof(OUT[i].y), :eL)
                sum_OUT = sum_OUT + OUT[i].y.eL
            elseif hasfield(typeof(OUT[i].y), :e)
                sum_OUT = sum_OUT + OUT[i].y.e
            end
            if haskey(OUT[i].Z,"eL")
                sum_OUT = sum_OUT + OUT[i].Z["eL"]
            end
            if haskey(OUT[i].Z,"e")
                sum_OUT = sum_OUT + OUT[i].Z["e"]
            end
        end
    end
    sum_e = sum_IN - sum_OUT
end
