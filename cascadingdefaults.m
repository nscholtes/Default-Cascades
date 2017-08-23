function [banks,CascadeOutput,Balancesheet_EV,banks_CDM_results_av] = ...
    cascadingdefaults(banks,n_banks,ActiveBanks,FailedBanks,adjacency_matrix,probability_matrix,ratios,shockpars,shockspecs)

%--------------------------------------------------------------------------
%% Initialisation and preallocation
%--------------------------------------------------------------------------

theta = ratios(1);
gamma = ratios(2);

shockmagnitude  = shockpars(1);
numshockedbanks = shockpars(2);

shocktype       = shockspecs{1};
targetcriterion = shockspecs{2};

%lossvec = zeros(1,n_banks);

init_assets  = zeros(1,n_banks);
init_capital = zeros(1,n_banks);
init_IBM_B   = zeros(1,n_banks);
init_IBM_L   = zeros(1,n_banks);

dispoutput = 0; % Show output on propagation dynamics (for debugging and future modifications)

IB_loan_matrix      = zeros(n_banks,n_banks);
IB_shockprop_matrix = zeros(n_banks,n_banks);

bank_capital = zeros(n_banks,1);

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

% (VI) Determining bank equity and deposits

for i = ActiveBanks
    banks(i).liabilities.total    = banks(i).assets.total;
    banks(i).liabilities.capital  = gamma*banks(i).liabilities.total;
    banks(i).liabilities.deposits = banks(i).assets.total - (banks(i).liabilities.capital+banks(i).liabilities.IB_Tot_Borrowing);   
end

% Collecting all borrower, lender and intermediary information

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

% Collecting initial values for normalisation of results from cascading defaults model
for i = ActiveBanks
    init_assets(i)  = banks(i).assets.total;
    init_capital(i) = banks(i).liabilities.capital;  
    init_IBM_B(i)   = banks(i).liabilities.IB_Tot_Borrowing;
end

Tot_init_assets   = sum(init_assets);
Tot_init_capital  = sum(init_capital);
Norm_init_capital = Tot_init_capital/Tot_init_assets;

Tot_init_IB_exp  = sum(init_IBM_B);
Norm_init_IB_exp = Tot_init_IB_exp/Tot_init_assets;

%--------------------------------------------------------------------------
%% Collecting individual bank information
%--------------------------------------------------------------------------

for i = ActiveBanks
    bank_capital(i) = banks(i).liabilities.capital;
end

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
    shockedbanks = randsample();
    n_runs      = 100;    
end

%shockedbanks = 1; %initial_shockedbank = shockedbanks;

TotFailedBanks      = zeros(1,n_runs);
Norm_TotFailedBanks = zeros(1,n_runs);

TotCapitalLoss      = zeros(1,n_runs);
Norm_TotCapitalLoss = zeros(1,n_runs);

shockterm = zeros(1,n_runs);

banks_CDM_results   = zeros(n_banks,3,n_runs);

for n = 1:n_runs   
    shockedbanks = shock_matrix(n,:);
    IB_loan_matrix      = zeros(n_banks,n_banks);
    IB_shockprop_matrix = zeros(n_banks,n_banks);
    
    [TotFailedBanks(n),Norm_TotFailedBanks(n),TotCapitalLoss(n),Norm_TotCapitalLoss(n),banks_CDM_results(:,:,n),shockterm(n)] =...
        cascadealgorithm(banks,n_banks,shockedbanks,ActiveBanks,FailedBanks,IB_loan_matrix,IB_shockprop_matrix,Tot_init_capital,shockpars,dispoutput);
end

% Compute averages across different runs of the cascade algorithm (with random shocks) across the same network

banks_CDM_results_av = mean(banks_CDM_results,3);

TFB_av       = mean(TotFailedBanks);
Norm_TFB_av  = mean(Norm_TotFailedBanks);
TCL_av       = mean(TotCapitalLoss);

Norm_TCL_av  = mean(Norm_TotCapitalLoss);
shockterm_av = mean(shockterm);

% Collecting output
CascadeOutput   = [TFB_av, Norm_TFB_av,TCL_av,Norm_TCL_av,shockterm_av];
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