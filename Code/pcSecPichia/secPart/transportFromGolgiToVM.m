function reaction = transportFromGolgiToVM(peptide,peptide_org)

%ALP pathway
reaction{1}.rxns = sprintf('%s_ALPtransport_sec_Apl6p_Aps3p_Apm3p_Apl5p_Vam3p_Clc1p_Chc1p_Arf1p_Swa2p_Vps1p_complex',peptide_org);


reaction{1}.rxnNames = sprintf('%s_Direct vacuol transit pathway_Apl6p_Aps3p_Apm3p_Apl5p_Vam3p_Clc1p_Chc1p_Arf1p_Swa2p_Vps1p',peptide);


reaction{1}.eq = sprintf('%s[g] + 4 gtp[c] + 4 h2o[c] => %s_folding[vm] + 4 gdp[c] + 4 pi[c] + 4 h[c]',peptide,peptide_org);

end