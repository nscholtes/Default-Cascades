function [banks,CascadeOutput,Balancesheet_EV,banks_CDM_results_av,banks_BS_results,shockedbanks_initBS_av,shockedbanks_centrality_av,initshock_magnitude_av] = ...
    cascadingdefaults(banks,n_banks,ActiveBanks,FailedBanks,adjacency_matrix,probability_matrix,...
    localnetworkmeasures,ratios,shockpars,shocktype,targetgroup,liquidityrisk,dispoutput)

%--------------------------------------------------------------------------
%% Initialisation and preallocation
%--------------------------------------------------------------------------

theta = ratios(1);
gamma = ratios(2);

tstop = 1000;

numshockedbanks = shockpars(2);

init_assets   =   zeros(n_banks,1);
init_ext_assets = zeros(n_banks,1);
init_deposits =   zeros(n_banks,1);
init_capital  =   zeros(n_banks,1);
init_IBM_B    =   zeros(n_banks,1);
init_IBM_L    =   zeros(n_banks,1);

IB_loan_matrix      = zeros(n_banks,n_banks);
IB_shockprop_matrix = zeros(n_banks,n_banks);

%--------------------------------------------------------------------------
%% Calibrating heterogenous bank balance sheets
%--------------------------------------------------------------------------

for i = ActiveBanks    
    
    banks(i).status = 1;
    
    % (I) Counterparty information from adjacency matrix 
    borrower_vec = adjacency_matrix(i,:);
    lender_vec   = adjacency_matrix(:,i);
    
    banks(i).borrower_ids = find(borrower_vec);
    banks(i).lender_ids   = find(lender_vec);
    
    banks(i).num_borrowers = numel(banks(i).borrower_ids);
    banks(i).num_lenders   = numel(banks(i).lender_ids);
    
    % (II) Balance sheet composition from bank size distribution and calibrated parameters
    
    banks(i).assets.externalassets = theta*banks(i).assets.total;
    banks(i).assets.IB_Tot_Lending = (1-theta)*banks(i).assets.total;
    
    banks(i).borrower_peculiarity = probability_matrix(i,banks(i).borrower_ids)';
    
    % (III) Interbank lending from (I) & (II)

    if banks(i).num_borrowers ~= 0
        banks(i).assets.IB_Bil_Lending = (banks(i).assets.IB_Tot_Lending*banks(i).borrower_peculiarity)./(sum(banks(i).borrower_peculiarity));
        
        IB_loan_matrix(i,banks(i).borrower_ids) = banks(i).assets.IB_Bil_Lending;
        banks(i).IB_IsLender = 1;
    else
        banks(i).assets.IB_Tot_Lending = 0;
        banks(i).assets.IB_Bil_Lending = 0;
        banks(i).IB_IsLender    = 0;
        
        banks(i).assets.externalassets = banks(i).assets.total;
    end        
end

%  (IV) Bank borrower information based on (III)

for i = ActiveBanks    
    if banks(i).num_lenders ~=0
        banks(i).liabilities.IB_Bil_Borrowing = IB_loan_matrix(banks(i).lender_ids,i);
        banks(i).liabilities.IB_Tot_Borrowing = sum(banks(i).liabilities.IB_Bil_Borrowing);
        banks(i).IB_IsBorrower = 1;
    else
        banks(i).liabilities.IB_Bil_Borrowing = 0;
        banks(i).liabilities.IB_Tot_Borrowing = 0;
        banks(i).IB_IsBorrower = 0;
    end   
    
    banks(i).IB_IsIntermediary = 0;
end

% (V) Classifying banks who borrow AND lend as intermediaries

for i = ActiveBanks
    if (banks(i).IB_IsLender == 1 && banks(i).IB_IsBorrower == 1)
        banks(i).IB_IsIntermediary = 1;
    end
end

% (VI) Determining bank capital and deposits

for i = ActiveBanks
    banks(i).liabilities.total    = banks(i).assets.total;
    banks(i).liabilities.capital  = gamma*banks(i).liabilities.total;
    banks(i).liabilities.deposits = banks(i).assets.total - (banks(i).liabilities.capital+banks(i).liabilities.IB_Tot_Borrowing);   
end

% (VII) Collecting all borrower, lender and intermediary information

for i = ActiveBanks    
    Borrower_vec(i)     = banks(i).IB_IsBorrower;
    Lender_vec(i)       = banks(i).IB_IsLender;
    Intermediary_vec(i) = banks(i).IB_IsIntermediary;    
end

All_Borrowers      = find(Borrower_vec);
All_Lenders        = find(Lender_vec);
All_Intermediaries = find(Intermediary_vec);

% Reorganising bank structure for clarity

permutationvector  = [1 9 3 7 5 4 6 8 10 11 2];
tempstore = banks;

clearvars banks
banks = orderfields(tempstore,permutationvector);

% (VIII) Collecting initial values
for i = ActiveBanks
    init_assets(i)   = banks(i).assets.total;
    init_deposits(i) = banks(i).liabilities.deposits;
    init_capital(i)  = banks(i).liabilities.capital;  
    init_IBM_B(i)    = banks(i).liabilities.IB_Tot_Borrowing;
    init_IBM_L(i)    = banks(i).assets.IB_Tot_Lending;
    init_ext_assets(i) = banks(i).assets.externalassets;
end

banks_BS_results = [init_assets init_capital init_IBM_B init_IBM_L];

Tot_init_assets     = sum(init_assets);
Tot_init_ext_assets = sum(init_ext_assets);

Tot_init_deposits   = sum(init_deposits);

Tot_init_capital  = sum(init_capital);
Norm_init_capital = Tot_init_capital/Tot_init_assets;

Tot_init_IB_exp  = sum(init_IBM_B);
Norm_init_IB_exp = Tot_init_IB_exp/Tot_init_assets;

banks_init = banks;

%--------------------------------------------------------------------------
%% Propagation dynamics
%--------------------------------------------------------------------------

%load('test.mat')

% (I) Initial shock to external assets

% Select set of shocked banks dependant on specification of random or targeted shock

if strcmp(shocktype,'random')
    n_runs       = 100; 
    shock_matrix = zeros(n_runs,numshockedbanks);
    for i = 1:n_runs
        shock_matrix(i,:) = randsample(All_Borrowers,numshockedbanks); % Select initially shocked bank(s) from amongst 
    end                                                                % all borrowers banks with at least one incoming edge
elseif strcmp(shocktype,'targeted')
    n_runs = 10;
    for i = 1:n_runs
        shock_matrix(i,:) = randsample(ActiveBanks(1:targetgroup),numshockedbanks);
    end  
end

% Matrix preallocation
TotFailedBanks       = zeros(1,n_runs);
Norm_TotFailedBanks  = zeros(1,n_runs);
TotCapitalLoss       = zeros(1,n_runs);
Norm_TotCapitalLoss  = zeros(1,n_runs);
TotDepositLoss      = zeros(1,n_runs);
Norm_TotDepositLoss = zeros(1,n_runs);

shockterm           = zeros(1,n_runs);
initshock_magnitude = zeros(1,n_runs);
banks_CDM_results   = zeros(n_banks,3,n_runs);
assetprice          = zeros(n_runs,tstop);


% Collect centralities (various measures) of initially shocked banks
shockedbanks_IDC = zeros(n_runs,numshockedbanks); shockedbanks_ODC = zeros(n_runs,numshockedbanks);
shockedbanks_IC  = zeros(n_runs,numshockedbanks); shockedbanks_OC  = zeros(n_runs,numshockedbanks);
shockedbanks_BC  = zeros(n_runs,numshockedbanks); shockedbanks_PRC = zeros(n_runs,numshockedbanks);

% Collect balance sheet properties of initially shocked banks    
shockedbanks_capital    = zeros(n_runs,numshockedbanks);
shockedbanks_assets     = zeros(n_runs,numshockedbanks);
shockedbanks_ext_assets = zeros(n_runs,numshockedbanks);
shockedbanks_IBL        = zeros(n_runs,numshockedbanks);
shockedbanks_IBB        = zeros(n_runs,numshockedbanks);

for n = 1:n_runs
    shockedbanks = shock_matrix(n,:);
    
    % Network measures
    shockedbanks_IDC(n,:) = localnetworkmeasures(shockedbanks,1)'; % In-degree centrality of shocked banks
    shockedbanks_ODC(n,:) = localnetworkmeasures(shockedbanks,2)'; % Out-degree centrality
    
    shockedbanks_IC(n,:)  = localnetworkmeasures(shockedbanks,3)'; % In-closeness centrality 
    shockedbanks_OC(n,:)  = localnetworkmeasures(shockedbanks,4)'; % Out-closeness centrality 
    
    shockedbanks_BC(n,:)  = localnetworkmeasures(shockedbanks,5)'; % Betweenness centrality
    shockedbanks_PRC(n,:) = localnetworkmeasures(shockedbanks,6)'; % PageRank centrality
    
    % Initial capital of shocked banks
    shockedbanks_capital(n,:)  = init_capital(shockedbanks);
    
    shockedbanks_assets(n,:)     = init_assets(shockedbanks);
    shockedbanks_ext_assets(n,:) = init_ext_assets(shockedbanks);
    
    shockedbanks_IBL(n,:)    = init_IBM_L(shockedbanks);
    shockedbanks_IBB(n,:)    = init_IBM_B(shockedbanks); 
    
    %IB_loan_matrix      = zeros(n_banks,n_banks);
    IB_shockprop_matrix = zeros(n_banks,n_banks);
    
    [banks,TotFailedBanks(n),Norm_TotFailedBanks(n),TotCapitalLoss(n),Norm_TotCapitalLoss(n),TotDepositLoss(n),Norm_TotDepositLoss(n),...
        banks_CDM_results(:,:,n),shockterm(n),initshock_magnitude(n),assetprice(n,1:shockterm(n)+1)] =...
        cascadealgorithm(banks_init,n_banks,shockedbanks,ActiveBanks,FailedBanks,probability_matrix,IB_loan_matrix,IB_shockprop_matrix,...
        Tot_init_capital,Tot_init_ext_assets,Tot_init_deposits,shockpars,dispoutput,liquidityrisk,tstop);
    
    d_assetprice_abs(n) = assetprice(n,1) - assetprice(n,shockterm(n)+1);                          % Absolute change in asset price
    d_assetprice_pct(n) =  ((assetprice(n,shockterm(n)+1) - assetprice(n,1))/assetprice(n,1))*100; % Percentage change in asset price
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute averages across different runs of the cascade algorithm (with random shocks) across the same network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Various centrality measures of shocked banks
shockedbanks_centrality = [shockedbanks_IDC shockedbanks_ODC shockedbanks_IC shockedbanks_OC shockedbanks_BC shockedbanks_PRC];

if numshockedbanks == 1
    shockedbanks_centrality_av = mean(shockedbanks_centrality);
else
    for i=1:size(localnetworkmeasures,2)
        SBC_temp(:,i) = mean(shockedbanks_centrality(:,(numshockedbanks*i)-(numshockedbanks-1):numshockedbanks*i),2);
    end
    shockedbanks_centrality_av = mean(SBC_temp);
end

shockedbanks_capital_av    = mean(mean(shockedbanks_capital,2));
shockedbanks_assets_av     = mean(mean(shockedbanks_assets,2));
shockedbanks_LR_av         = mean(mean(shockedbanks_capital./shockedbanks_assets,2)); % Shocked banks' Leverage ratio = capital/assets
shockedbanks_ext_assets_av = mean(mean(shockedbanks_ext_assets,2));
shockedbanks_IBB_av        = mean(mean(shockedbanks_IBB,2)); 
shockedbanks_IBL_av        = mean(mean(shockedbanks_IBL,2));

shockedbanks_initBS_av = [shockedbanks_capital_av shockedbanks_assets_av shockedbanks_LR_av shockedbanks_ext_assets_av...
    shockedbanks_IBB_av shockedbanks_IBL_av];

initshock_magnitude_av = mean(initshock_magnitude); % Magnitude of initial shock

banks_CDM_results_av = mean(banks_CDM_results,3);   % Cascading default model results for individual banks (bank level)

% Cascading default model results over all banks (network level)
TFB_av       = mean(TotFailedBanks);
Norm_TFB_av  = mean(Norm_TotFailedBanks);

TCL_av       = mean(TotCapitalLoss);
Norm_TCL_av  = mean(Norm_TotCapitalLoss);

TDL_av       = mean(TotDepositLoss);
Norm_TDL_av  = mean(Norm_TotDepositLoss);

shockterm_av = mean(shockterm);

% Change in asset price over the simulation
d_assetprice_abs_av = mean(d_assetprice_abs);
d_assetprice_pct_av = mean(d_assetprice_pct);
d_assetprice_av     = [d_assetprice_abs_av d_assetprice_pct_av];

% Collecting output at the network level
CascadeOutput   = [TFB_av,Norm_TFB_av,TCL_av,Norm_TCL_av,TDL_av,Norm_TDL_av,shockterm_av d_assetprice_av];
Balancesheet_EV = [Tot_init_assets, Tot_init_capital, Norm_init_capital Tot_init_IB_exp Norm_init_IB_exp];

%--------------------------------------------------------------------------
% Plot cascading default dynamics
%--------------------------------------------------------------------------
% 
%  t  = 0:0.1:10;
%  y1 = sin(t);
%  y2 = sin(t + 0.1) + 0.5;
%  y3 = cos(t).^2;
%  y4 = cos(t + 0.2).^2 - 0.2;
% 
%  figure
%  hold on
%  fill([t fliplr(t)],[y2 fliplr(y1)] ,'r')
%  fill([t fliplr(t)],[y4 fliplr(y3)] ,'b')
%  grid
%  alpha(0.15)

end