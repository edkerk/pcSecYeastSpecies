%% setMedia
function model = setMediaPP(model,type)

% type  = 1: minimal media (Delft media) (default)
%       = 2: yeast nitrogen base without amino acids
%       = 3: YNB+CSM-Ura (DOI: 10.1039/B803529F)
%       = 4: SD_2_SCAA
%       = 5: yeast nitrogen base with all amino acids

exchangeRxns = findExcRxns(model);
model.lb(exchangeRxns) = 0;
model.ub(exchangeRxns) = 1000;


% 1: minimal media (Delft media)                
desiredExchanges_1 = {'Ex_nh4'; ... % ammonium exchange
                    'Ex_o2'; ... % oxygen exchange
                    'Ex_pi'; ... % phosphate exchange
                    'Ex_so4'; ... % sulphate exchange
                    'Ex_fe2'; ... % iron exchange, for test of expanded biomass def
                    'Ex_h'; ... % hydrogen exchange
                    'Ex_h2o'; ... % water exchange
                    'Ex_na1'; ... % sodium exchange
                    'Ex_k'};    % potassium exchange
% 2: yeast nitrogen base without amino acids
desiredExchanges_2 = {'Ex_btn'; ... % biotin exchange(C10H16N2O3S)
                      'Ex_thm'; ... % thiamine exchange(C12H17N4OS+)
                      'Ex_nac'; ... % nicotinate exchange(C6H4NO2-)
                      'Ex_inost'; ... % myo-inositol exchange(C6H12O6)
                      'Ex_4abz'; ... % 4-aminobenzoate exchange(C7H6NO2-)
                      'Ex_pnto_R'; ... % (R)-pantothenate exchange(C9H17NO5)
                      'Ex_chol'; ... % 'Choline exchange       
                      'Ex_ribflv'};    % riboflavin exchange(C17H20N4O6)
% 3: CSM-Ura
CSM_Ura_Exchanges =  {'Ex_his_L'; ... % L-histidine exchange
                      'Ex_met_L'; ... % L-methionine exchange
                      'Ex_trp_L'; ... % L-tryptophan exchange
                      'Ex_ade'; ... % adenine exchange
                      'Ex_arg_L'; ... % L-arginine exchange
                      'Ex_asp_L'; ... % L-aspartate exchange
                      'Ex_ile_L'; ... % L-isoleucine exchange
                      'Ex_leu_L'; ... % L-leucine exchange
                      'Ex_lys_L'; ... % L-lysine exchange
                      'Ex_phe_L'; ... % L-phenylalanine exchange
                      'Ex_thr_L'; ... % L-threonine exchange
                      'Ex_tyr_L'; ... % L-tyrosine exchange
                      'Ex_val_L'};    % L-valine exchange

% 4: minimal + AA(14)
SD_2_SCAA             =    {'Ex_arg_L'; ...% L-arginine exchange
                            'Ex_asp_L';...%L-aspartate exchange
                            'Ex_glu_L';...%L-glutamate exchange
                            'Ex_gly';...%L-glycine exchange
                            'Ex_his_L';...%L-histidine exchange
                            'Ex_ile_L';...%L-isoleucine exchange
                            'Ex_leu_L';...%L-leucine exchange
                            'Ex_lys_L';...%L-lysine exchange
                            'Ex_met_L';...%L-methionine exchange
                            'Ex_phe_L';...%L-phenylalanine exchange
                            'Ex_thr_L';...%L-threonine exchange
                            'Ex_trp_L';...%L-tryptophan exchange
                            'Ex_tyr_L';...%L-tyrosine exchange
                            'Ex_val_L';...%L-valine exchange
                            'Ex_ura'};  %uracil exchange
% 5 YNB with aa
YNB                    =    {'Ex_arg_L'; ...% L-arginine exchange
                            'Ex_asp_L';...%L-aspartate exchange
                            'Ex_glu_L';...%L-glutamate exchange
                            'Ex_gly';...%L-glycine exchange
                            'Ex_his_L';...%L-histidine exchange
                            'Ex_ile_L';...%L-isoleucine exchange
                            'Ex_leu_L';...%L-leucine exchange
                            'Ex_lys_L';...%L-lysine exchange
                            'Ex_met_L';...%L-methionine exchange
                            'Ex_phe_L';...%L-phenylalanine exchange
                            'Ex_thr_L';...%L-threonine exchange
                            'Ex_trp_L';...%L-tryptophan exchange
                            'Ex_tyr_L';...%L-tyrosine exchange
                            'Ex_val_L';...%L-valine exchange
                            'Ex_ala_L';...%L-alanine exchange
                            'Ex_asp_L';...%L-asparagine exchange
                            'Ex_cys_L';...%L-cysteine exchange
                            'Ex_gln_L';...%L-glutamine exchange
                            'Ex_pro_L';...%L-proline exchange
                            'Ex_ser_L' }; %L-serine exchange
% 7 SC selection medium
SC                    =    {'Ex_arg_L';...% L-arginine exchange
                            'Ex_asp_L';...%L-aspartate exchange
                            'Ex_glu_L';...%L-glutamate exchange
                            'Ex_his_L';...%L-histidine exchange
                            'Ex_ile_L';...%L-isoleucine exchange
                            'Ex_leu_L';...%L-leucine exchange
                            'Ex_lys_L';...%L-lysine exchange
                            'Ex_met_L';...%L-methionine exchange
                            'Ex_phe_L';...%L-phenylalanine exchange
                            'Ex_ser_L';...%L-serine exchange
                            'Ex_thr_L';...%L-threonine exchange
                            'Ex_trp_L';...%L-tryptophan exchange
                            'Ex_tyr_L';...%L-tyrosine exchange
                            'Ex_val_L';...%L-valine exchange
                            'Ex_ura'};  %uracil exchange

% 8 YEP medium
YEP                    =    {'Ex_arg_L'; ...% L-arginine exchange
                            'Ex_asn_L';...%L-asparagine exchange
                            'Ex_glu_L';...%L-glutamate exchange
                            'Ex_gly';...%L-glycine exchange
                            'Ex_his_L';...%L-histidine exchange
                            'Ex_ile_L';...%L-isoleucine exchange
                            'Ex_leu_L';...%L-leucine exchange
                            'Ex_lys_L';...%L-lysine exchange
                            'Ex_met_L';...%L-methionine exchange
                            'Ex_phe_L';...%L-phenylalanine exchange
                            'Ex_thr_L';...%L-threonine exchange
                            'Ex_trp_L';...%L-tryptophan exchange
                            'Ex_tyr_L';...%L-tyrosine exchange
                            'Ex_val_L';...%L-valine exchange
                            'Ex_ala_L';...%L-alanine exchange
                            'Ex_asp_L';...%L-asparagine exchange
                            'Ex_cys_L';...%L-cysteine exchange
                            'Ex_gln_L';...%L-glutamine exchange
                            'Ex_pro_L';...%L-proline exchange
                            'Ex_ser_L';...%L-serine exchange
                            'Ex_thymd';...%thymidine exchange
                            'Ex_ura'};%Uracil exchange

% % %                             'Ex_dad_2';...%2''-deoxyadenosine exchange
% % %                             'Ex_dgsn';...%2''-deoxyguanosine exchange
% % %                             'Ex_thymd';...%thymidine exchange
% % %                             'Ex_dcyt'};%deoxycytidine exchange

          
                           


uptakeRxnIndexes     = findRxnIDs(model,desiredExchanges_1);
uptakeRxnIndexes_2     = findRxnIDs(model,desiredExchanges_2);

% % blockedRxnIndex      = findRxnIDs(model,blockedExchanges);
CSM_Ura_RxnIndex     = findRxnIDs(model,CSM_Ura_Exchanges);
SD_2_SCAA_RxnIndex = findRxnIDs(model,SD_2_SCAA);
YNB_RxnIndex = findRxnIDs(model,YNB);
SC_RxnIndex = findRxnIDs(model,SC);
YEP_RxnIndex = findRxnIDs(model,YEP);



if length(find(uptakeRxnIndexes~= 0)) ~= length(desiredExchanges_1)
    warning('Not all exchange reactions were found.');
end

model.lb(uptakeRxnIndexes(uptakeRxnIndexes~=0)) = -1000;
model.ub(uptakeRxnIndexes(uptakeRxnIndexes~=0)) = 1000;

% % % model.lb(blockedRxnIndex) = 0;
% % % model.ub(blockedRxnIndex) = 0;

if type == 2
     model.lb(uptakeRxnIndexes_2) = -2; 
elseif type == 3
     model.lb(uptakeRxnIndexes_2) = -2; 
     model.lb(CSM_Ura_RxnIndex) = -0.08; % max coefficiency of AA in biomass times max mu
elseif type == 4
     model.lb(uptakeRxnIndexes_2) = -2; 
     model.lb(SD_2_SCAA_RxnIndex) = -0.08; % max coefficiency of AA in biomass times max mu
elseif type == 5
     model.lb(uptakeRxnIndexes_2) = -2; 
     model.lb(YNB_RxnIndex) = -0.1; % max coefficiency of AA in biomass times max mu
elseif type == 7
     model.lb(uptakeRxnIndexes_2) = -2; 
     model.lb(SC_RxnIndex) = -0.08; % max coefficiency of AA in biomass times max mu
elseif type == 8
     model.lb(uptakeRxnIndexes_2) = -2; 
     model.lb(YEP_RxnIndex) = -0.08; % max coefficiency of AA in biomass times max mu
    

end

