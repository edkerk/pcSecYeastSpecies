function reaction = transportFromGolgiToS(peptide,peptide_org)

%HDSV pathway
reaction{1}.rxns = sprintf('%s_HDSVI_sec_Arf1p_Pep12p_Swa2p_Chc1p_Clc1p_Apl4p_Apl2p_Apm1p_Aps1p_complex',peptide_org);
reaction{2}.rxns = sprintf('%s_HDSVII_sec_Vps1p_Chc1p_Clc1p_complex',peptide_org);


reaction{1}.rxnNames = sprintf('%s_HDSVI_Arf1p_Pep12p_Swa2p_Chc1p_Clc1p_Apl4p_Apl2p_Apm1p_Aps1p',peptide);
reaction{2}.rxnNames = sprintf('%s_HDSVII_Vps1p_Chc1p_Clc1p',peptide);


reaction{1}.eq = sprintf('%s[g] + gtp[c] + h2o[c] => %s[ce] + gdp[c] + pi[c] + h[c]',peptide,peptide);
reaction{2}.eq = sprintf('%s[ce] + gtp[c] + h2o[c] => %s_folding[e] + gdp[c] + pi[c] + h[c]',peptide,peptide_org);

end