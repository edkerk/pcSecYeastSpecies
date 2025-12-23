function NGAM_T = getNGAMT(t)
    % T is in K, a single value
    T = t + 273.15;
    % Nested function to calculate NGAM based on temperature T
    function result = NGAM_function(T)
        result = 0.740 + 5.893/(1+exp(31.920-(T-273.15))) + 6.12e-6*(T-273.15-16.72)^4;
    end

    % Define lower and upper bounds in Kelvin
    lb = 10 + 273.15;
    ub = 60 + 273.15;
    
    % Calculate NGAM_T based on the temperature T
    if T < lb
        NGAM_T = NGAM_function(lb);
    elseif T > ub
        NGAM_T = NGAM_function(ub);
    else
        NGAM_T = NGAM_function(T);
    end
end




