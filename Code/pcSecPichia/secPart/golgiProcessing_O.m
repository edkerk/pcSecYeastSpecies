function [newModel,peptide_name,rxns] = golgiProcessing_O(model,peptide_h,peptide_org,Length_total,OG,onlyrxns)
rxns =[];
if OG >0
    Length = Length_total/40;

    reaction{1}.rxns = sprintf('%s_GLOG_Golgi_O_linked_manosylation_I_sec_KTR_complex',peptide_org);
    
    reaction{1}.rxnNames = sprintf('%s_GLOG_Golgi_O_linked_manosylation_I_sec_KTR_complex',peptide_h);
    
    reaction{1}.eq = sprintf('%s[g] + %d gdpmann[g] => %s_GOG_G1[g] + %d gdp[g]',peptide_h,2*OG,peptide_h,2*OG);%%%according to Radoman et al,2020

        if onlyrxns == 1
            rxns = [rxns;{reaction{1}.rxns}];
        else
            model=addYeastReaction(model,reaction{1}.eq,{reaction{1}.rxns},{reaction{1}.rxnNames});
            rxns = [rxns;{reaction{1}.rxns}];
        end
    newModel = model;
    peptide_name = sprintf('%s_GOG_G1',peptide_h);
else
    newModel = model;
    peptide_name = peptide_h;
    
end