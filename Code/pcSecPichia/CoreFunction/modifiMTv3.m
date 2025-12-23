%Rxns Added to generate new model

function model= modifiMTv3(model)
%GTP/GDP & gthox&gthrd in ER
[model,rxnIDexists] = addReaction(model,{'GTPter','GTP transport into ER'} ,'gtp[c] <=> gtp[er]');
[model,rxnIDexists] = addReaction(model,{'GDPter','GDP transport into ER'} ,'gdp[c] <=> gdp[er]');
[model,rxnIDexists] = addReaction(model,{'GTHOXter','gthox transport into ER'} ,'gthox[c] <=> gthox[er]');
[model,rxnIDexists] = addReaction(model,{'GTHRDter','gthrd transport into ER'} ,'gthrd[c] <=> gthrd[er]');

[model,rxnIDexists] = addReaction(model,{'ERAD4M','Ubiquitin of protein'} ,'h2o[c] + atp[c] + Ubiquitin[c] -> h[c] + ppi[c] + amp[c] + Ubiquitin_for_Transfer[c]');
[model,rxnIDexists] = addReaction(model,{'PDI_reoxidation_GSSG','PDI_reoxidation_GSSG'} ,'PDI-ox[er] + 2 gthrd[er] -> PDI[er] + gthox[er]');
[model,rxnIDexists] = addReaction(model,{'PDI_reoxidation_ERO1','PDI_reoxidation_ERO1'} ,'PDI[er] + o2[er] -> PDI-ox[er] + h2o2[er]');

%Modify the model

[model,rxnIDexists] = addReaction(model,{'hexccoater','hexccoa transport into ER'} ,'hexccoa[c] <=> hexccoa[er]');%%%a problem in original model

%GPI synthesis
[model,rxnIDexists] = addReaction(model,{'ethampter','Ethanolamine phosphate transport into ER'} ,'ethamp[c] <=> ethamp[er]');
[model,rxnIDexists] = addReaction(model,{'hdcater','Hexadecanoate transport into ER'} ,'hdca[c] <=> hdca[er]');
[model,rxnIDexists] = addReaction(model,{'oleateter','oleate transport into ER'} ,'ocdcea[c] <=> ocdcea[er]');
[model,rxnIDexists] = addReaction(model,{'cer3_26ter','Ceramide-3 (Phytosphingosine:n-C26:0OH) transport into ER'} ,'cer3_26[c] <=> cer3_26[er]');
[model,rxnIDexists] = addReaction(model,{'N-Acetyle-Glucoseamine-transferase','N-Acetyle-Glucoseamine-transferase'} ,'ptd1ino_SC[c] + uacgam[c] -> h[c] + udp[c] + 6-(N-acetyl-alpha-D-glucosaminyl)-1-phosphatidyl-1D-myo-inositol[c]');
[model,rxnIDexists] = addReaction(model,{'GPI-anchor assembly, step 2','GPI-anchor assembly, step 2'} ,'6-(N-acetyl-alpha-D-glucosaminyl)-1-phosphatidyl-1D-myo-inositol[c] + h2o[c] -> 6-(alpha-D-glucosaminyl)-1-phosphatidyl-1D-myo-inositol[c] + ac[c] + h[c]');
[model,rxnIDexists] = addReaction(model,{'GPI-anchor assembly','GPI-anchor assembly'} ,'6-(alpha-D-glucosaminyl)-1-phosphatidyl-1D-myo-inositol[c] -> 6-(alpha-D-glucosaminyl)-1-phosphatidyl-1D-myo-inositol[er]');
[model,rxnIDexists] = addReaction(model,{'acylation of GPI inositol at 2 position, GPI-anchor assembly, step 3','acylation of GPI inositol at 2 position, GPI-anchor assembly, step 3'} ,'6-(alpha-D-glucosaminyl)-1-phosphatidyl-1D-myo-inositol[er] + pmtcoa[er] -> 6-(alpha-D-glucosaminyl)-O-acyl-1-phosphatidyl-1D-myo-inositol[er] + coa[er]');
[model,rxnIDexists] = addReaction(model,{'GPI-anchor assembly, step 5','GPI-anchor assembly, step 5'} ,'6-(alpha-D-glucosaminyl)-O-acyl-1-phosphatidyl-1D-myo-inositol[er] + dolmanp[er] -> 6-O-(alpha-D-mannosyl-(1-4)-alpha-D-glucosaminyl)-O-acyl-1-phosphatidyl-1D-myo-inositol[er] + dolp[er]');
[model,rxnIDexists] = addReaction(model,{'GPI-anchor assembly, step 6','GPI-anchor assembly, step 6'} ,'6-O-(alpha-D-mannosyl-(1-4)-alpha-D-glucosaminyl)-O-acyl-1-phosphatidyl-1D-myo-inositol[er] + pe_SC[er] -> 6-O-{2-O-((2-aminoethyl)phosphoryl)-alpha-D-mannosyl-(1-4)-alpha-D-glucosaminyl}-O-acyl-1-phosphatidyl-1D-myo-inositol[er] + 12dgr_SC[er]');
[model,rxnIDexists] = addReaction(model,{'GPI-anchor assembly, step 7','GPI-anchor assembly, step 7'} ,'6-O-{2-O-((2-aminoethyl)phosphoryl)-alpha-D-mannosyl-(1-4)-alpha-D-glucosaminyl}-O-acyl-1-phosphatidyl-1D-myo-inositol[er] + dolmanp[er] -> 6-O-{alpha-D-mannosyl-(1-6)-2-O-((2-aminoethyl)phosphoryl)-alpha-D-mannosyl-(1-4)-alpha-D-glucosaminyl}-O-acyl-1-phosphatidyl-1D-myo-inositol[er] + dolp[er]');
[model,rxnIDexists] = addReaction(model,{'GPI-anchor assembly, step 8','GPI-anchor assembly, step 8'} ,'6-O-{alpha-D-mannosyl-(1-6)-2-O-((2-aminoethyl)phosphoryl)-alpha-D-mannosyl-(1-4)-alpha-D-glucosaminyl}-O-acyl-1-phosphatidyl-1D-myo-inositol[er] + dolmanp[er] -> 6-O-{alpha-D-mannosyl-(1-2)-alpha-D-mannosyl-(1-6)-2-O-((2-aminoethyl)phosphoryl)-alpha-D-mannosyl-(1-4)-alpha-D-glucosaminyl}-O-acyl-1-phosphatidyl-1D-myo-inositol[er] + dolp[er]');
[model,rxnIDexists] = addReaction(model,{'GPI-anchor assembly, step 9','GPI-anchor assembly, step 9'} ,'6-O-{alpha-D-mannosyl-(1-2)-alpha-D-mannosyl-(1-6)-2-O-((2-aminoethyl)phosphoryl)-alpha-D-mannosyl-(1-4)-alpha-D-glucosaminyl}-O-acyl-1-phosphatidyl-1D-myo-inositol[er] + dolmanp[er] -> 6-O-alpha-D-mannosyl-(1-2)-{alpha-D-mannosyl-(1-2)-alpha-D-mannosyl-(1-6)-2-O-((2-aminoethyl)phosphoryl)-alpha-D-mannosyl-(1-4)-alpha-D-glucosaminyl}-O-acyl-1-phosphatidyl-1D-myo-inositol[er] + dolp[er]');
[model,rxnIDexists] = addReaction(model,{'GPI-anchor assembly, step 10','GPI-anchor assembly, step 10'} ,'6-O-alpha-D-mannosyl-(1-2)-{alpha-D-mannosyl-(1-2)-alpha-D-mannosyl-(1-6)-2-O-((2-aminoethyl)phosphoryl)-alpha-D-mannosyl-(1-4)-alpha-D-glucosaminyl}-O-acyl-1-phosphatidyl-1D-myo-inositol[er] + pe_SC[er] -> 6-O-alpha-D-mannosyl-(1-2)-{alpha-D-mannosyl-2-O-((2-aminoethyl)phosphoryl)-(1-2)-alpha-D-mannosyl-(1-6)-2-O-((2-aminoethyl)phosphoryl)-alpha-D-mannosyl-(1-4)-alpha-D-glucosaminyl}-O-acyl-1-phosphatidyl-1D-myo-inositol[er] + 12dgr_SC[er]');
[model,rxnIDexists] = addReaction(model,{'GPI-anchor assembly, step 11','GPI-anchor assembly, step 11'} ,'6-O-alpha-D-mannosyl-(1-2)-{alpha-D-mannosyl-2-O-((2-aminoethyl)phosphoryl)-(1-2)-alpha-D-mannosyl-(1-6)-2-O-((2-aminoethyl)phosphoryl)-alpha-D-mannosyl-(1-4)-alpha-D-glucosaminyl}-O-acyl-1-phosphatidyl-1D-myo-inositol[er] + pe_SC[er] -> 6-O-2-O-((2-aminoethyl)phosphoryl)-alpha-D-mannosyl-(1-2)-{alpha-D-mannosyl-2-O-((2-aminoethyl)phosphoryl)-(1-2)-alpha-D-mannosyl-(1-6)-2-O-((2-aminoethyl)phosphoryl)-alpha-D-mannosyl-(1-4)-alpha-D-glucosaminyl}-O-acyl-1-phosphatidyl-1D-myo-inositol[er] + 12dgr_SC[er]');


