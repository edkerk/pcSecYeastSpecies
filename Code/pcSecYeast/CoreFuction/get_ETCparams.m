% % get_ETCparams
function ETCparams = get_ETCparams(Tm_info,t)
    % Calculate the fraction of enzyme in native state
    % dHTH, dSTS are enthalpy and entropy at TH and TS, respectively
    % dCpu is the heat capacity change upon unfolding
    % T is the temperature, can be a single value or an array, in K
    
    R = 8.314;  % J/(mol*K)
    TH = 373.5;  % Reference temperature in Kelvin
    TS = 385;    % Another reference temperature in Kelvi
    T = t + 273.15;
    ETCparams = Tm_info;


    for i = 1:length(ETCparams.Gene)
        Tm = ETCparams.Tm(i)+273.15;
        use_length_based_algorithm = true; 

        if isfield(ETCparams, 'T90') && ~isnan(ETCparams.T90(i)) && ETCparams.T90(i) > ETCparams.Tm(i)
            T90 = ETCparams.T90(i) + 273.15; 
            slope = 299.58;
            intercept = 20008;
    
            % Define the matrix A and vector B for the system of equations
            A = [1, -slope, 0;
                1, -Tm, Tm - TH - Tm * log(Tm / TS);
                1, -T90, T90 - TH - T90 * log(T90 / TS)];
     
            B = [intercept; 0; -R * T90 * log(9)];
    
         % Solve the system of equations
            solution = A \ B;
    
            % Extract the results
            dHTH = solution(1);
            dSTS = solution(2);
            dCpu = solution(3);
% %             if dHTH >= 0 && dSTS >= 0 && dCpu >= 0
            if dCpu >= 0
                use_length_based_algorithm = false;
            end
        end
        if use_length_based_algorithm
            N = strlength(ETCparams.Sequence(i));
        % Define the gas constant
                % Define the function to calculate deltaG(Tm)
    
            dHTH =  (4*N + 143) * 1000;
            dSTS = 13.27 * N + 448;
            
            func = @(dCp) dHTH + dCp * (Tm - TH) - Tm * dSTS - Tm * dCp * log(Tm / TS);
            
            % Use fsolve to solve for dCpu

            options = optimset('Display', 'off');  % Turn off fsolve output
            dCpu = fsolve(func, 10000, options);
        end

        % Calculate the Gibbs free energy change (deltaG) at temperature T
        dGu = dHTH + dCpu * (T - TH) - T * dSTS - T * dCpu * log(T / TS);
        
        % Calculate the fraction of enzyme in native state
        fNT = 1 ./ (1 + exp(-dGu ./ (R .* T)));
        fUT = 1-fNT;


        
        ETCparams.dCpu(i,1) = dCpu;
        ETCparams.dHTH(i,1) = dHTH;
        ETCparams.dSTS(i,1) = dSTS;
        ETCparams.dGu(i,1) = dGu;
        ETCparams.fNT(i,1) = fNT;
        ETCparams.fUT(i,1) = fUT;
    end
end


%% get_dH_dS_dCpu_params

function [dHTH,dSTS,dCpu] = get_dH_dS_dCpu_params(Tm, T90) 
    % With knowing Tm and T90, get dHTH, dSTS and dCpu by solving
    % dHTH = slope*dSTS+intercept
    % deltaG(Tm) = 0
    % deltaG(T90) = -RTln9
    % 
    % to get dSTS, dCpu.
    % 
    % Tm, T90 are in K
    % dHTH, is in J/mol
    % dSTS is in J/mol/K


    TH = 373.5;
    TS = 385;
    R = 8.314;
    slope = 299.58;
    intercept = 20008;
    
    % Define the matrix A and vector B
    A = [1, -slope, 0;
         1, -Tm, Tm - TH - Tm*log(Tm/TS);
         1, -T90, T90 - TH - T90*log(T90/TS)];
     
    B = [intercept; 0; -R*T90*log(9)];
    
    % Solve the system of linear equations
    X = A \ B;
    
    % Assign the results
    dHTH = X(1);
    dSTS = X(2);
    dCpu = X(3);
end

%% get_dH_dS_dCpu_from_TmLength

function [dHTH, dSTS, dCpu] = get_dH_dS_dCpu_from_TmLength(Tm, N)
    % In case of negative obtained from get_dH_dS_dCpu_from_TmT90(Tm, T90), or 
    % if there is no T90 data available, use the same equations from Sawle 
    % and Ghosh, Biophysical Journal, 2011 for deltaH* and deltaS*.
    % Then calculate dCpu by solving deltaG(Tm) = 0.
    %
    % Tm is in K
    % dHTH is in J/mol
    % dSTS is in J/mol/K
    
    % Calculate dHTH and dSTS using the given formulas
    dHTH = (4*N + 143) * 1000;
    dSTS = 13.27 * N + 448;
    
    % Define constants
    TH = 373.5;  % Reference temperature (K)
    TS = 385;    % Another reference temperature (K)
    
    % Define the function to calculate deltaG(Tm)
    func = @(dCp) dHTH + dCp * (Tm - TH) - Tm * dSTS - Tm * dCp * log(Tm / TS);
    
    % Use fsolve to solve for dCpu
    dCpu = fsolve(func, 10000);
    
%     % Return the calculated values
%     return
end

%% get_dGu

function dGu = get_dGu(T, dHTH, dSTS, dCpu)
    % calculate the deltaG of unfolding process at temperature T
    % dHTH, dSTS are enthalpy and entropy at TH and TS, respectively
    % dCpu is the heat capacity change upon unfolding
    
    % Define the reference temperatures
    TH = 373.5;  % Reference temperature in Kelvin
    TS = 385;    % Another reference temperature in Kelvin

    % Calculate the Gibbs free energy change (deltaG) at temperature T
    dGu = dHTH + dCpu * (T - TH) - T * dSTS - T * dCpu * log(T / TS);
end

%% get_fNT

function f = get_fNT(T, dHTH, dSTS, dCpu)
    % Calculate the fraction of enzyme in native state
    % dHTH, dSTS are enthalpy and entropy at TH and TS, respectively
    % dCpu is the heat capacity change upon unfolding
    % T is the temperature, can be a single value or an array, in K
    
    % Define the gas constant
    R = 8.314;  % J/(mol*K)
    
    % Calculate the Gibbs free energy change (dGu) at temperature T
    dGu = get_dGu(T, dHTH, dSTS, dCpu);
    
    % Calculate the fraction of enzyme in native state
    f = 1 ./ (1 + exp(-dGu ./ (R .* T)));
end


