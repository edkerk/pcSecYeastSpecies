function [ETCparamsT,enzymedataT] = calculate_kcatT(ETCparams, enzymedata, strain,t)
    % Using Transition state theory to calculate kcat at temperature T.
    % dHTH, dSTS: entropy and enthalpy at convergence temperatures. Protein
    % unfolding process.
    % dCpu, heat capacity change upon unfolding.
    % kcatTopt: kcat values at optimal temperature
    % Topt, optimal temperature of the enzyme, in K
    % T, temperature, in K

    % Constants
    R = 8.314;
    TH = 373.5;
    TS = 385;
    T0 = 30 + 273.15;
    T = t + 273.15;
    dCpt = -6300;
    
    if strcmp(strain, 'KM')
        Topt = 38 + 273.15;
    elseif  strcmp(strain, 'PP')
        Topt = 28 + 273.15;
     else  
        Topt = 30 + 273.15;
     end
    % Initialize empty cell array 
    all_genes = {};
    all_kcat = [];
    ETCparams.kcatTopt = zeros(length(ETCparams.Gene),1);
    % Iterate over each row in the enzymedata.subunit matrix
    for row = 1:size(enzymedata.subunit, 1)
        current_genes = enzymedata.subunit(row, ~cellfun('isempty', enzymedata.subunit(row, :)));
        all_genes = [all_genes; current_genes'];
        all_kcat = [all_kcat; repmat(enzymedata.kcat(row), length(current_genes), 1)];
    end
    
    % Iterate over each gene in the all_genes cell array
    for k = 1:length(all_genes)
        gene = all_genes{k};
        kcat = all_kcat(k);
      
        % Find the indices in ETCparams where the gene matches
        match_idx = find(strcmp(ETCparams.Gene, gene));
        
        % If there are matching indices
        if ~isempty(match_idx)
            % Update the kcatTopt column in ETCparams for each matching index
            for j = 1:length(match_idx)
                ETCparams.kcatTopt(match_idx(j),1) = kcat;
            end
        end
        
    end
    nonzero_kcat = ETCparams.kcatTopt(ETCparams.kcatTopt ~= 0);
    mean_kcat = mean(nonzero_kcat);
    ETCparams.kcatTopt(ETCparams.kcatTopt == 0) = mean_kcat;

    for i = 1:length(ETCparams.kcatTopt)
        dCpu = ETCparams.dCpu(i,1);
        dHTH = ETCparams.dHTH(i,1);
        dSTS = ETCparams.dSTS(i,1);
        kcatTopt = ETCparams.kcatTopt(i,1);
        % Use the equation from solvedHT.m and re-organized
        dGuTopt = dHTH + dCpu * (Topt - TH) - Topt * dSTS - Topt * dCpu * log(Topt / TS);
        dHt = dHTH + dCpu * (Topt - TH) - dCpt * (Topt - T0) - R * Topt - ...
              (dHTH + dCpu * (Topt - TH)) / (1 + exp(-dGuTopt / (R * Topt)));
    
        % Calculate kcat at reference Temperature
        kcat0 = kcatTopt / exp(log(Topt / T0) - (dHt + dCpt * (Topt - T0)) / (R * Topt) + ...
               dHt / (R * T0) + dCpt * log(Topt / T0) / R);
    
        % Calculate kcat at given temperature
        kcatT = kcat0 * exp(log(T / T0) - (dHt + dCpt * (T - T0)) / (R * T) + ...
                dHt / (R * T0) + dCpt * log(T / T0) / R);
        
        ETCparams.kcat0(i,1) = kcat0;
        ETCparams.kcatT(i,1) = kcatT;
    end
    genenumber = min(length(ETCparams.Gene),length(enzymedata.subunit));
    for n = 1:genenumber
        gene = ETCparams.Gene(n);
        match_idx = find(strcmp(enzymedata.subunit(:,1), gene));
        if ~isempty(match_idx)          
            for m = 1:length(match_idx)
                enzymedata.kcat(match_idx(m),1) = ETCparams.kcatT(n);
            end
        end

    end
    enzymedata.kcat(enzymedata.kcat == 0 | isnan(enzymedata.kcat)) = mean_kcat;


    enzymedataT = enzymedata;
    ETCparamsT = ETCparams;

end


% ETCparams = table(ETCparams.Gene,ETCparams.Sequence,ETCparams.Tm,ETCparams.dCpu,ETCparams.dHTH,ETCparams.dSTS,ETCparams.dGu,ETCparams.fNT,ETCparams.fUT,ETCparams.kcatTopt,ETCparams.kcat0,ETCparams.kcatT);
% ETCparams.Properties.VariableNames = { ...
% 'Gene', 'Sequence', 'Tm', 'dCpu', 'dHTH', ...
% 'dSTS', 'dGu', 'fNT', 'fUT', 'kcatTopt', ...
% 'kcat0', 'kcatT'};
% writetable(ETCparams, 'ETCparams.csv');