function [newModel,peptide_name,rxns] = coat_trans_membrane(model,peptide,peptide_org,GPI,trans,onlyrxns)
rxns = [];
if trans >0 && GPI == 0 

    reaction{1}.rxns = sprintf('%s_COPII_TransM_ERGL1B_sec_Sec12p_Sar1p_Sec23p_Sec24p_Bet1p_Bos1p_complex',peptide_org);
    reaction{2}.rxns = sprintf('%s_COPII_TransM_ERGL1B_sec_Sec12p_Sar1p_Shl23p_Lst1p_Bet1p_Bos1p_complex',peptide_org);
    reaction{3}.rxns = sprintf('%s_COPII_ERGL_sec_Sec13p_Sec31p_Sec16p_Sed4p_Sec5p_Sec17p_complex',peptide_org);
    reaction{4}.rxns = sprintf('%s_COPII_ERGL_sec_Ypt1p_Uso1p_bug1p_Bet3p_Bet5p_Trs20p_Trs23p_Trs31p_Trs33p_complex',peptide_org);
    
    reaction{1}.rxnNames = sprintf('%s_COPII_TransM_ERGL1B_sec_Sec12p_Sar1p_Sec23p_Sec24p_Bet1p_Bos1p_complex Pre budding complex forming for soluble proteins',peptide);
    reaction{2}.rxnNames = sprintf('%s_COPII_TransM_ERGL1B_sec_Sec12p_Sar1p_Shl23p_Lst1p_Bet1p_Bos1p_complex Pre budding complex forming for soluble proteins',peptide);
    reaction{3}.rxnNames =  sprintf('%s_COPII_common_ERGLA_Sec12p_Sar1p_Sec23p_Sec24p_Erv29p COPII formation',peptide);
    reaction{4}.rxnNames = sprintf('%s_COPII_common_ERGLA_Ypt1p_Uso1p_bug1p_Bet3p_Bet5p_Trs20p_Trs23p_Trs31p_Trs33p COPII fusion',peptide);
      
    reaction{1}.eq = sprintf('%s[er] + gtp[er] + h2o[er] => %s_COP_coated[er] + gdp[er] + h[er] + pi[er]',peptide,peptide);
    reaction{2}.eq = sprintf('%s[er] + gtp[er] + h2o[er] => %s_COP_coated[er] + gdp[er] + h[er] + pi[er]',peptide,peptide);
    reaction{3}.eq = sprintf('%s_COP_coated[er] => %s_COP_coated[c]',peptide,peptide);
    reaction{4}.eq = sprintf('%s_COP_coated[c] + gtp[c] + h2o[c] => %s[g] + gdp[c] + h[c] + pi[c]',peptide,peptide);
    

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