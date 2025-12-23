function reaction = Post_Translation_translocation(peptide,Length)
%This function makes five reactions that translocates a peptide using post
%translation
%The function returns these five reaction

reaction{1}.rxns = sprintf('%s_Post_translation_PSTA_sec_RAC_complex',peptide);
reaction{2}.rxns = sprintf('%s_Post_translation_PSTA_sec_Ssa1_Ydj1_Snl1_complex',peptide);
reaction{3}.rxns = sprintf('%s_Post_translation_PSTA_sec_SEC61SEC63C_complex',peptide);
reaction{4}.rxns = sprintf('%s_Post_translation_PSTA_sec_BIP_NEFS_complex',peptide);
reaction{5}.rxns = sprintf('%s_Post_translation_TC_sec_SPC_complex',peptide);
reaction{6}.rxns = sprintf('%s_export_sp_to_c',peptide);

reaction{1}.rxnNames = sprintf('%s_Post_translation_PSTA_sec_RAC_complex',peptide);
reaction{2}.rxnNames = sprintf('%s_Post_translation_PSTA_sec_Ssa1_Ydj1_Snl1_complex',peptide);
reaction{3}.rxnNames = sprintf('%s_Post_translation_PSTA_sec_SEC61SEC63C_complex',peptide);
reaction{4}.rxnNames = sprintf('%s_Post_translation_PSTA_sec_BIP_NEFS_complex',peptide);
reaction{5}.rxnNames = sprintf('%s_Signal peptidase',peptide);
reaction{6}.rxnNames = sprintf('%s_export sp to cytosol',peptide);

reaction{1}.eq = sprintf('%s_peptide[c] => %s_translocate_1[c]',peptide,peptide);
reaction{2}.eq= sprintf('%s_translocate_1[c] + atp[c] + h2o[c] => %s_translocate_2[c] + adp[c] + h[c] + pi[c]',peptide,peptide);
reaction{3}.eq= sprintf('%s_translocate_2[c] => %s_translocate_3[c]',peptide,peptide);
reaction{4}.eq= sprintf('%s_translocate_3[c] + %.15f atp[c] + %.15f h2o[c] => %s[er] + %.15f adp[c] + %.15f h[c] + %d pi[c]',peptide,Length,Length,peptide,Length,Length,Length);
reaction{5}.eq= sprintf('%s_translocate_3[c] + h2o[c] => %s[er] + %s_sp[er]',peptide,peptide,peptide);
reaction{6}.eq = sprintf('%s_sp[er] => %s_sp[c]',peptide,peptide);
end


