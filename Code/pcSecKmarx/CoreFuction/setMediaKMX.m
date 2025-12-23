%% setMedia
function model = setMediaKMX(model,type)

% type  = 1: minimal KMXmedia (Delft media) (default)
%       = 2: yeast nitrogen base without amino acids
%       = 3: YNB+CSM-Ura (DOI: 10.1039/B803529F)
%       = 4: SD_2_SCAA
%       = 5: yeast nitrogen base with all amino acids
%       = 6：KMX iSM996 test media 
%       = 7: SC selection medium


exchangeRxns = findExcRxns(model);
model.lb(exchangeRxns) = 0;
model.ub(exchangeRxns) = 1000;

blockedExchanges = {'r_1877'...% pyridoxal exchange
                    'r_1771'} ; % bicarbonate exchange
% 1: minimal media (Delft media)                
desiredExchanges_1 = {'r_1727'; ... % ammonium exchange
                    'r_1725'; ... % oxygen exchange
                    'r_1729'; ... % phosphate exchange
                    'r_1728'; ... % sulphate exchange
                    'r_1731'; ... % iron exchange, for test of expanded biomass def
                    'r_1724'; ... % hydrogen exchange
                    'r_1723'; ... % water exchange
                    'r_1884'; ... % sodium exchange
                    'r_1875'};    % potassium exchange
% 2: yeast nitrogen base without amino acids
desiredExchanges_2 = {'r_1772'; ... % biotin exchange(C10H16N2O3S)
                      'r_1734'; ... % thiamine exchange(C12H17N4OS+)
                      'r_1733'; ... % pyridoxine exchange >high uptake<(C8H11NO3)
                      'r_1732'; ... % nicotinate exchange(C6H4NO2-)
                      'r_1730'; ... % myo-inositol exchange(C6H12O6)
                      'r_1736'; ... % 4-aminobenzoate exchange(C7H6NO2-)
                      'r_1735'; ... % (R)-pantothenate exchange(C9H17NO5)
                      'r_1754'; ... %% folic acid exchange(5-formyltetrahydrofolic acid exchange)(C20H23N7O7)
                      'r_1881'};    % riboflavin exchange(C17H20N4O6)
% 3: CSM-Ura
CSM_Ura_Exchanges =  {'r_1840'; ... % L-histidine exchange
                      'r_1845'; ... % L-methionine exchange
                      'r_1853'; ... % L-tryptophan exchange
                      'r_1760'; ... % adenine exchange
                      'r_1831'; ... % L-arginine exchange
                      'r_1833'; ... % L-aspartate exchange
                      'r_1842'; ... % L-isoleucine exchange
                      'r_1843'; ... % L-leucine exchange
                      'r_1844'; ... % L-lysine exchange
                      'r_1847'; ... % L-phenylalanine exchange
                      'r_1852'; ... % L-threonine exchange
                      'r_1854'; ... % L-tyrosine exchange
                      'r_1855'};    % L-valine exchange

% 4: minimal + AA(14)
SD_2_SCAA             =    {'r_1831'; ...% L-arginine exchange
                            'r_1833';...%L-aspartate exchange
                            'r_1837';...%L-glutamate exchange
                            'r_1839';...%L-glycine exchange
                            'r_1840';...%L-histidine exchange
                            'r_1842';...%L-isoleucine exchange
                            'r_1843';...%L-leucine exchange
                            'r_1844';...%L-lysine exchange
                            'r_1845';...%L-methionine exchange
                            'r_1847';...%L-phenylalanine exchange
                            'r_1852';...%L-threonine exchange
                            'r_1853';...%L-tryptophan exchange
                            'r_1854';...%L-tyrosine exchange
                            'r_1855';...%L-valine exchange
                            'r_1897'};  %uracil exchange
% 5 YNB with aa
YNB                    =    {'r_1831'; ...% L-arginine exchange
                            'r_1833';...%L-aspartate exchange
                            'r_1837';...%L-glutamate exchange
                            'r_1839';...%L-glycine exchange
                            'r_1840';...%L-histidine exchange
                            'r_1842';...%L-isoleucine exchange
                            'r_1843';...%L-leucine exchange
                            'r_1844';...%L-lysine exchange
                            'r_1845';...%L-methionine exchange
                            'r_1847';...%L-phenylalanine exchange
                            'r_1852';...%L-threonine exchange
                            'r_1853';...%L-tryptophan exchange
                            'r_1854';...%L-tyrosine exchange
                            'r_1855';...%L-valine exchange
                            'r_1828';...%L-alanine exchange
                            'r_1832';...%L-asparagine exchange
                            'r_1835';...%L-cysteine exchange
                            'r_1838';...%L-glutamine exchange
                            'r_1848';...%L-proline exchange
                            'r_1850' }; %L-serine exchange
%6 KMX iSM996 test media 
KM                     =    {'r_1723';...%H2O exchange
                             'r_1724';...%H+ exchange
                             'r_1725';...%oxygen exchange
                             'r_1726';...%D-glucose exchange
                             'r_1727';...%ammonium exchange
                             'r_1728';...%sulphate exchange       
                             'r_1729';...%phosphate exchange
                             'r_1730';...%myo-inositol exchange
                             'r_1731';...%iron exchange
                             'r_1732';...%nicotinate exchange
                             'r_1733';...%pyridoxine exchange
                             'r_1734';...%thiamine exchange
                             'r_1735';...%(R)-pantothenate exchange
                             'r_1736'};%4-aminobenzoate exchange

% 7 SC selection medium
SC                    =    {'r_1831';...% L-arginine exchange
                            'r_1833';...%L-aspartate exchange
                            'r_1837';...%L-glutamate exchange
                            'r_1840';...%L-histidine exchange
                            'r_1842';...%L-isoleucine exchange
                            'r_1843';...%L-leucine exchange
                            'r_1844';...%L-lysine exchange
                            'r_1845';...%L-methionine exchange
                            'r_1847';...%L-phenylalanine exchange
                            'r_1850';...%L-serine exchange
                            'r_1852';...%L-threonine exchange
                            'r_1853';...%L-tryptophan exchange
                            'r_1854';...%L-tyrosine exchange
                            'r_1855';...%L-valine exchange
                            'r_1897'};  %uracil exchange




uptakeRxnIndexes     = findRxnIDs(model,desiredExchanges_1);
uptakeRxnIndexes_2     = findRxnIDs(model,desiredExchanges_2);

blockedRxnIndex      = findRxnIDs(model,blockedExchanges);
CSM_Ura_RxnIndex     = findRxnIDs(model,CSM_Ura_Exchanges);
SD_2_SCAA_RxnIndex = findRxnIDs(model,SD_2_SCAA);
YNB_RxnIndex = findRxnIDs(model,YNB);
KM_RxnIndex = findRxnIDs(model,KM);
SC_RxnIndex = findRxnIDs(model,SC);


if length(find(uptakeRxnIndexes~= 0)) ~= length(desiredExchanges_1)
    warning('Not all exchange reactions were found.');
end

model.lb(uptakeRxnIndexes(uptakeRxnIndexes~=0)) = -1000;

model.lb(blockedRxnIndex) = 0;
model.ub(blockedRxnIndex) = 0;


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
     model.lb(YNB_RxnIndex) = -0.08; % max coefficiency of AA in biomass times max mu
elseif type == 6
     model.lb(KM_RxnIndex(KM_RxnIndex~=0)) = -1000;
     model.ub(KM_RxnIndex(KM_RxnIndex~=0)) = 1000;
elseif type == 7
     model.lb(uptakeRxnIndexes_2) = -2; 
     model.lb(SC_RxnIndex) = -0.08; % max coefficiency of AA in biomass times max mu
end
