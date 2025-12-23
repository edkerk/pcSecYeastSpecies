function  [newModel,peptide_name,rxns] = coat_GPI(model,peptide,peptide_org,GPI,onlyrxns)
rxns = [];
if GPI > 0
    reaction{1}.rxns = sprintf('%s_COPII_GPI_ERGL1C_sec_Sec12p_Sar1p_Sec23p_Sec24p_Emp24p_Erp1p_Erp2p_Erv25p_Bos1p_Bet1p_complex',peptide_org);
    reaction{2}.rxns = sprintf('%s_COPII_GPI_ERGL1C_sec_Sec12p_Sar1p_Shl23p_Lst1p_Emp24p_Erp1p_Erp2p_Erv25p_Bos1p_Bet1p_complex',peptide_org);%Sec23/24 homologous Shl23/Lst1 was added
    reaction{3}.rxns = sprintf('%s_COPII_ERGL_sec_Sec13p_Sec31p_Sec16p_Sed4p_Sec5p_Sec17p_complex',peptide_org);
    reaction{4}.rxns = sprintf('%s_COPII_ERGL_sec_Ypt1p_Uso1p_Bet3p_Bet5p_Trs20p_Trs23p_Trs31p_Trs33p_complex',peptide_org);
    
    reaction{1}.rxnNames = sprintf('%s_COPII_GPI_ERGL1C_sec_Sec12p_Sar1p_Sec23p_Sec24p_Emp24p_Erp1p_Erp2p_Erv25p_Bos1p_Bet1p_complex Pre budding complex forming for GPI-anchored proteins',peptide);
    reaction{2}.rxnNames = sprintf('%s_COPII_GPI_ERGL1C_sec_Sec12p_Sar1p_Shl23p_Lst1p_Emp24p_Erp1p_Erp2p_Erv25p_Bos1p_Bet1p_complex Pre budding complex forming for GPI-anchored proteins',peptide);
    reaction{3}.rxnNames = sprintf('%s_COPII_GPI_ERGLB_Sec13p_Sec31p_Sec16p_Sed4p_Sec5p_Sec17p COPII formation',peptide);
    reaction{4}.rxnNames = sprintf('%s_COPII_GPI_ERGLB_Ypt1p_Uso1p_Bet3p_Bet5p_Trs20p_Trs23p_Trs31p_Trs33p_complex COPII fusion',peptide); 
    reaction{1}.eq = sprintf('%s[er] + gtp[er] + h2o[er] => %s_GPI_G1_COP[er] + gdp[er] + h[er] + pi[er]',peptide,peptide);
    reaction{2}.eq = sprintf('%s[er] + gtp[er] + h2o[er] => %s_GPI_G1_COP[er] + gdp[er] + h[er] + pi[er]',peptide,peptide);
    reaction{3}.eq = sprintf('%s_GPI_G1_COP[er] => %s_GPI_G2_COP[c]',peptide,peptide);
    reaction{4}.eq = sprintf('%s_GPI_G2_COP[c] + gtp[c] + h2o[c] => %s[g] + gdp[c] + h[c] + pi[c]',peptide,peptide);
    
    for i=1:4
        if onlyrxns == 1
            rxns = [rxns;{reaction{i}.rxns}];
        else
            model=addYeastReaction(model,reaction{i}.eq,{reaction{i}.rxns},{reaction{i}.rxnNames});
            rxns = [rxns;{reaction{i}.rxns}];
        end
    end
    newModel = model;
    peptide_name = sprintf('%s',peptide);
else
    newModel = model;
    peptide_name = peptide;
    
end