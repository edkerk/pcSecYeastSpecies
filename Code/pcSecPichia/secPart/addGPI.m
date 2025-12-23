function [newModel peptide_name,rxns] = addGPI(model,peptide,peptide_org,Length_total,GPI,onlyrxns)
rxns = [];
if GPI >0
    if GPI > 1
        GPI = 1;
    end
    reaction{1}.rxns = sprintf('%s_GPIRI_sec_GPIR_complex',peptide_org);
    reaction{2}.rxns = sprintf('%s_GPIRII_sec_Bst1p_complex',peptide_org);
    reaction{3}.rxns = sprintf('%s_GPIRIII_sec_Per1p_complex',peptide_org);
    reaction{4}.rxns = sprintf('%s_GPIRIV_sec_Gup1p_complex',peptide_org);
    reaction{5}.rxns = sprintf('%s_GPIRV_sec_Cwh43p_Las21p_Mcd4p_complex',peptide_org);
    reaction{6}.rxns = sprintf('%s_GPIRVI_sec_Ted1p_complex',peptide_org);
    reaction{7}.rxns = sprintf('%s_GPIRIB_sec_GPIR_complex',peptide_org);

    reaction{1}.rxnNames = sprintf('%s_GPI_the luminal part of the protein is attached to a GPI anchor via a phosphoethanolaminethe luminal part of the protein is attached to a GPI anchor via a phosphoethanolamine_GPIR_complex',peptide);
    reaction{2}.rxnNames = sprintf('%s_GPI_removes the acyl chain from the inositol_Bst1p',peptide);
    reaction{3}.rxnNames = sprintf('%s_GPI_removal of the unsaturated acyl chain at the sn-2 position of diacylglycerol to form lyso-GPI_Per1p',peptide);
    reaction{4}.rxnNames = sprintf('%s_GPI_C26 saturated acyl chain is transferred to the sn-2 position_Gup1p',peptide);
    reaction{5}.rxnNames = sprintf('%s_GPI_the lipid moiety is further changed to ceramide consisting of PHS with a hydroxy-C26 fatty acid_Cwh43p_Las21p_Mcd4p_complex',peptide);
    reaction{6}.rxnNames = sprintf('%s_GPI_removes a PEtN on the second mannose_Ted1p',peptide);
    reaction{7}.rxnNames = sprintf('%s_GPI_GPI transfer using sugar and lipid dispatched from misfolded GPI protein_GPIR_complex',peptide);

    reaction{1}.eq = sprintf('%s[er] + %.15f 6-O-2-O-((2-aminoethyl)phosphoryl)-alpha-D-mannosyl-(1-2)-{alpha-D-mannosyl-2-O-((2-aminoethyl)phosphoryl)-(1-2)-alpha-D-mannosyl-(1-6)-2-O-((2-aminoethyl)phosphoryl)-alpha-D-mannosyl-(1-4)-alpha-D-glucosaminyl}-O-acyl-1-phosphatidyl-1D-myo-inositol[er] => %s_GPI_G1[er] + %.15f h2o[er]',peptide,GPI,peptide,GPI);
    reaction{2}.eq = sprintf('%s_GPI_G1[er] => %s_GPI_G2[er] + %.15f hdca[er]',peptide,peptide,GPI);
    reaction{3}.eq = sprintf('%s_GPI_G2[er] + %.15f h2o[er] => %s_GPI_G3[er] + %.15f ocdcea[er]',peptide,GPI,peptide,GPI); 
    reaction{4}.eq = sprintf('%s_GPI_G3[er] + %.15f hexccoa[er] => %s_GPI_G4[er] + %.15f coa[er] + %.15f h[er]',peptide,GPI,peptide,GPI,GPI);
    reaction{5}.eq = sprintf('%s_GPI_G4[er] + %.15f cer3_26[er] => %s_GPI_G5[er] + %.15f 12dgr_SC[er]',peptide,GPI,peptide,GPI);
    reaction{6}.eq = sprintf('%s_GPI_G5[er] + %.15f h2o[er] => %.15f ethamp[er] + %s_GPI_G6[er]',peptide,GPI,GPI,peptide);
    %the last one is for using misfolded protein degradation product
    reaction{7}.eq = sprintf('%s[er] + %.15f 6-O-2-O-alpha-D-mannosyl-(1-2)-{alpha-D-mannosyl-2-O-((2-aminoethyl)phosphoryl)-(1-2)-alpha-D-mannosyl-(1-6)-2-O-((2-aminoethyl)phosphoryl)-alpha-D-mannosyl-(1-4)-alpha-D-glucosaminyl}-O-inositol-P-ceramide C (C26)[er] => %s_GPI_G6[er] + %.15f h2o[er]',peptide,GPI,peptide,GPI);



    for i=1:7
         if onlyrxns == 1
            rxns = [rxns;{reaction{i}.rxns}];
         else
            model=addYeastReaction(model,reaction{i}.eq,{reaction{i}.rxns},{reaction{i}.rxnNames});
            rxns = [rxns;{reaction{i}.rxns}];
         end
    end
    newModel = model;
    peptide_name = sprintf('%s_GPI_G6',peptide);
else
    newModel = model;
    peptide_name = peptide;

end
