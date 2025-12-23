function Tm_info = mapTmparams(Tm_pred,enzymedata,protein_info)
    % Initialize output structure
    Tm_info = struct();
    Tm_info.Gene = enzymedata.proteins;
    Tm_info.Tm = zeros(length(Tm_info.Gene),1);

    % Extract protein identifiers and sequences from protein_info
    protein = protein_info(2:end,2);
    sequence = protein_info(2:end,11);

    % Map sequences from protein_info to the corresponding genes in enzymedata
    for i = 1:length(protein)
        match_idx = find(strcmp(Tm_info.Gene, protein(i)));
        if ~isempty(match_idx)
            seq_clean = extractBefore(sequence(i), "*"); % remove stop codon marker (*)
            for j = 1:length(match_idx)
                Tm_info.Sequence(match_idx(j),1) = seq_clean;
            end
        end
    end

    % Map predicted Tm/Topt values from Tm_pred to corresponding genes
    for i = 1:length(Tm_info.Gene)
        Gene = Tm_info.Gene(i);
        match_idx = find(strcmp(Tm_pred.Gene,Gene));
        if ~isempty(match_idx)
            Tm_info.Tm(i,1)   = max(Tm_pred.predict_Tm(match_idx));
            Tm_info.Topt(i,1) = max(Tm_pred.Topt(match_idx));
        end
    end
end
