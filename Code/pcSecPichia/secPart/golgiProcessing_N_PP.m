function [newModel,peptide_name,rxns] = golgiProcessing_N_PP(model,peptide,peptide_org,Length_total,NG,onlyrxns)
rxns = [];
if NG >0
    Length = Length_total/40;

        reaction{1}.rxns = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_I_sec_Och1p_complex',peptide_org);%Both Man9 and hypermannosylation
        reaction{2}.rxns = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_II_sec_MPOLI_complex',peptide_org);%~20% of hypermannosylation
        reaction{3}.rxns = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_III_sec_MPoLII_complex',peptide_org);
        reaction{4}.rxns = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_II_sec_Mnn2pA_complex',peptide_org);
        reaction{5}.rxns = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_II_sec_Mnn2pB_complex',peptide_org);
        reaction{6}.rxns = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_II_sec_Mnn2pC_complex',peptide_org);%there are 3 homologous of Mnn2p in P.pastoris
    
        reaction{1}.rxnNames = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_I_sec_Och1p_complex',peptide);
        reaction{2}.rxnNames = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_II_sec_MPOLI_complex',peptide);
        reaction{3}.rxnNames = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_III_sec_MPoLII_complex',peptide);
        reaction{4}.rxnNames = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_II_sec_Mnn2pA_complex',peptide);
        reaction{5}.rxnNames = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_II_sec_Mnn2pB_complex',peptide);
        reaction{6}.rxnNames = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_II_sec_Mnn2pC_complex',peptide);
    
    
        reaction{1}.eq = sprintf('%s[g] + %d gdpmann[g] => %s_GNG_G1[g] + %d gdp[g]',peptide,NG,peptide,NG);
        reaction{2}.eq=  sprintf('%s_GNG_G1[g] + %d gdpmann[g] => %s_GNG_G2[g] + %d gdp[g]',peptide,9*0.2*NG,peptide,9*0.2*NG);
        reaction{3}.eq=  sprintf('%s_GNG_G2[g] + %d gdpmann[g] => %s_GNG_G3[g] + %d gdp[g]',peptide,30*0.2*NG,peptide,30*0.2*NG);
        reaction{4}.eq=  sprintf('%s_GNG_G3[g] + %d gdpmann[g] => %s_GNG_G4[g] + %d gdp[g]',peptide,0.2*NG,peptide,0.2*NG);
        reaction{5}.eq=  sprintf('%s_GNG_G3[g] + %d gdpmann[g] => %s_GNG_G4[g] + %d gdp[g]',peptide,0.2*NG,peptide,0.2*NG);
        reaction{6}.eq=  sprintf('%s_GNG_G3[g] + %d gdpmann[g] => %s_GNG_G4[g] + %d gdp[g]',peptide,0.2*NG,peptide,0.2*NG);

        for i=1:6
            if onlyrxns == 1
                rxns = [rxns;{reaction{i}.rxns}];
            else
                model=addYeastReaction(model,reaction{i}.eq,{reaction{i}.rxns},{reaction{i}.rxnNames});
                rxns = [rxns;{reaction{i}.rxns}];
            end
        end

    newModel = model;
    peptide_name = sprintf('%s_GNG_G4',peptide);
else
    newModel = model;
    peptide_name = peptide;
    
end        
    
% %     newModel = model;
% %     peptide_name_noh = sprintf('%s_GNG_G4',peptide);
% %     peptide_name_h = peptide;  %%%%%%%
% % else
% %     newModel = model;
% %     peptide_name_noh = peptide;
% %     peptide_name_h = peptide;
% %     
% % end