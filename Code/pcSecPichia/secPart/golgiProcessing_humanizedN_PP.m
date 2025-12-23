function [newModel,peptide_name,rxns] = golgiProcessing_humanizedN_PP(model,peptide,peptide_org,Length_total,NG,onlyrxns)
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
    

        reaction{1} .eq=  sprintf('%s[g] + %d h2o[g] => %s_M5[g] + %d man[g]',peptide,3*NG,peptide,3*NG);  %% Man I        
        reaction{2} .eq=  sprintf('%s_M5[g] + %d h2o[g] => %s_M3[g] + %d man[g]',peptide,NG,peptide,NG);  %% Man II
        reaction{3} .eq=  sprintf('%s_M3[g] + %d uacgam[g] => %s_GNG_M3GNA3[g] + %d h[g] + %d udp[g]',peptide,NG,peptide,NG,NG);  %% GlcNAcT-I(GnTI)
        reaction{4} .eq=  sprintf('%s_GNG_M3GNA3[g] + %d uacgam[g] => %s_GNG_M3GNA4[g] + %d h[g] + %d udp[g]',peptide,NG,peptide,NG,NG);  %% GlcNAcT-II(GnTII)
        reaction{5} .eq=  sprintf('%s_GNG_M3GNA4[g] + %d udpgal[g]  => %s_GNG_M3GNA4G2[g] + %d udp[g]',peptide,2*NG,peptide,2*NG);  %% GalT
        reaction{6} .eq=  sprintf('%s_GNG_M3GNA4G2[g] + %d cmpacnam[g] => %s_GNG_M3GNA4G2S2[g] + %d cmp[g]',peptide,2*NG,peptide,2*NG);  %% SiaT

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



