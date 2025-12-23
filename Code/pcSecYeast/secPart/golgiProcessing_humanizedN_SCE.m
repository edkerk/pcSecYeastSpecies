function [newModel,peptide_name,rxns] = golgiProcessing_humanizedN_SCE(model,peptide,peptide_org,Length_total,NG,onlyrxns)
rxns = [];
if NG >0
    Length = Length_total/40;

        reaction{1}.rxns = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_I_sec_ManI_complex',peptide_org);%Both Man9 and hypermannosylation
        reaction{2}.rxns = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_II_sec_ManII_complex',peptide_org);%~20% of hypermannosylation
        reaction{3}.rxns = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_III_sec_GnTI_complex',peptide_org);
        reaction{4}.rxns = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_II_sec_GnTII_complex',peptide_org);
        reaction{5}.rxns = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_II_sec_GalT_complex',peptide_org);
        reaction{6}.rxns = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_II_sec_SiaT_complex',peptide_org);%there are 3 homologous of Mnn2p in P.pastoris
    
        reaction{1}.rxnNames = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_I_sec_ManI_complex',peptide);
        reaction{2}.rxnNames = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_II_sec_ManII_complex',peptide);
        reaction{3}.rxnNames = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_III_sec_GnTI_complex',peptide);
        reaction{4}.rxnNames = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_II_sec_GnTII_complex',peptide);
        reaction{5}.rxnNames = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_II_sec_GalT_complex',peptide);
        reaction{6}.rxnNames = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_II_sec_SiaT_complex',peptide);
    

        reaction{1} .eq=  sprintf('%s[g] + %d H2O [Golgi] => %s_M5[g] + %d D-mannose [Golgi]',peptide,3*NG,peptide,3*NG);  %% Man I        
        reaction{2} .eq=  sprintf('%s_M5[g] + %d H2O [Golgi] => %s_M3[g] + %d D-mannose [Golgi]',peptide,NG,peptide,NG);  %% Man II
        reaction{3} .eq=  sprintf('%s_M3[g] + %d UDP-N-acetyl-alpha-D-glucosamine [Golgi] => %s_GNG_M3GNA3[g] + %d H+ [Golgi] + %d UDP [Golgi]',peptide,NG,peptide,NG,NG);  %% GlcNAcT-I(GnTI)
        reaction{4} .eq=  sprintf('%s_GNG_M3GNA3[g] + %d UDP-N-acetyl-alpha-D-glucosamine [Golgi] => %s_GNG_M3GNA4[g] + %d H+ [Golgi] + %d UDP [Golgi]',peptide,NG,peptide,NG,NG);  %% GlcNAcT-II(GnTII)
        reaction{5} .eq=  sprintf('%s_GNG_M3GNA4[g] + %d UDP-D-galactose [Golgi]  => %s_GNG_M3GNA4G2[g] + %d UDP [Golgi]',peptide,2*NG,peptide,2*NG);  %% GalT
        reaction{6} .eq=  sprintf('%s_GNG_M3GNA4G2[g] + %d CMP-N-acetylneuraminate [Golgi] => %s_GNG_M3GNA4G2S2[g] + %d CMP [Golgi]',peptide,2*NG,peptide,2*NG);  %% SiaT

% % 
        for i=1:6
            if onlyrxns == 1
                rxns = [rxns;{reaction{i}.rxns}];
            else
                model=addYeastReaction(model,reaction{i}.eq,{reaction{i}.rxns},{reaction{i}.rxnNames});
                rxns = [rxns;{reaction{i}.rxns}];
            end
        end
    
    newModel = model;
% % 
    peptide_name = sprintf('%s_GNG_M3GNA4G2S2',peptide);
else
    newModel = model;
    peptide_name = peptide;
    
end




% % 
% % idx = find(strcmp(model.metNames,'KLMA_40391_M8_M3[Golgi]'));
% % model.mets(idx)
% % Tag = findRxnsFromMets(model, ans)
% % printRxnFormula(model, 'metNameFlag',true, 'rxnAbbr', Tag);
% % printRxnFormula(model,Tag);
% % 
% % idx = find(strcmp(model.mets,'KLMA_40391_M8_M3[g]'));
% % model.metNames(idx)
% % 
% % Rxns = Mets2Rxns(model, 'KLMA_40391_M8_GNG_M3GNA4G2S2[g]');
% % Rxns.formulas
% % 
% % Rxns = Mets2Rxns(model1, 'PAS_chr2-2_0107_M8_GNG_G4[g]');
% % Rxns.formulas
% % printRxnFormula(NGmodel, 'metNameFlag',true, 'rxnAbbr', Rxns.rxns);
% % 
% % 
% % idx = find(strcmp(model.metNames,'G10694'));
% % nonZeroRowsIdx = find(full(model.S(idx,:)));
% % relatedReactions = model.rxns(nonZeroRowsIdx);
% % printRxnFormula(model, 'metNameFlag',true, 'rxnAbbr', relatedReactions);
% % 
% % 
% % printRxnFormula(NGmodel,



