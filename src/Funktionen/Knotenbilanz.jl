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