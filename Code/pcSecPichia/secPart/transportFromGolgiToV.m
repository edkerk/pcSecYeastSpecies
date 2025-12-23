function reaction = transportFromGolgiToV(peptide,peptide_org)

%CPY pathway
reaction{1}.rxns = sprintf('%s_CPYI_sec_Gga1p_Gga2p_Arf1p_Apl4p_Apl2p_Apm1p_Aps1p_Chc1p_Clc1p_Pep12p_Vps45p_Vps5p_Swa2p_complex',peptide_org);
reaction{2}.rxns = sprintf('%s_CPYII_sec_Vps4p_Vps27p_Apl6p_Aps3p_Apm3p_Apl5p_Vam3p_complex',peptide_org);


reaction{1}.rxnNames = sprintf('%s_CPYI_Gga1p_Gga2p_Arf1p_Apl4p_Apl2p_Apm1p_Aps1p_Chc1p_Clc1p_Pep12p_Vps45p_Vps5p_Swa2p',peptide);
reaction{2}.rxnNames = sprintf('%s_CPYII_Vps4p_Vps27p_Apl6p_Aps3p_Apm3p_Apl5p_Vam3p',peptide);


reaction{1}.eq = sprintf('%s[g] + 4 gtp[c] + 4 h2o[c] => %s_CPY_G1[v] + 4 gdp[c] + 4 pi[c] + 4 h[c]',peptide,peptide);
reaction{2}.eq=  sprintf('%s_CPY_G1[v] + atp[c] + h2o[c] => %s_folding[v] + adp[c] + h[c] + pi[c]',peptide,peptide_org);

end