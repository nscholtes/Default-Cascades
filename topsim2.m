function [LHD,Data_table,bank_results,Networks] = topsim2()

%--------------------------------------------------------------------------
% MATLAB code for 'Default cascades and systemic risk on different interbank network topologies'
% by Nicolas K. Scholtes (2017)

% Copyright (C) Nicolas K. Scholtes, 2017
% Distributed under GPL 3.0

addpath(genpath('/Users/nscholte/Desktop/Research/Ch.3 - Systemic risk/MATLAB Codes'));

% Output figures to directory with LaTeX draft for automatic updating
fig_output  = '/Users/nscholte/Desktop/Research/Ch.3 - Systemic risk/Drafts/Figures/';
% Results tables loaded into R for empirical treatment
data_output = '/Users/nscholte/Desktop/Research/Ch.3 - Systemic risk/R Codes/';
tab_output  = '/Users/nscholte/Desktop/Research/Ch.3 - Systemic risk/Drafts/Tables/';

%--------------------------------------------------------------------------
%% CALIBRATION OF SIMULATION PARMETERS
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bank size/fitness distribution parameters
%%% Featured in Equations () and () of the paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_banks = 250; % Number of banks

%Bank size distribution
a_min   = [10,50]; % Domain: Size of smallest bank
a_ratio = [20,50]; % Domain: Ratio of largest to smallest bank
gamma_a = [2 3];   % Power law exponent for bank size distribution     

% Disassortative network generation algorithm parameters
alpha   = [0.5 1];
beta    = [0.5 1];
d       = [0.1 1]; % Calibrate network density

networkpars = [a_min; a_ratio; gamma_a; alpha; beta; d];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial balance sheet weights
%%% Featured in Equations () and () of the paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta = [0.5,0.9];   % external asset/asset ratio
gamma = [0.01,0.1] ; % capital/asset ratio

BSpars = [theta; gamma];

%--------------------------------------------------------------------------
% Latin hypercube sampling
%--------------------------------------------------------------------------
numsamples = 1000;
LH_input   = [networkpars; BSpars];
[LHD]      = createLH(LH_input,numsamples);
%--------------------------------------------------------------------------
% Generating the network
%--------------------------------------------------------------------------

tic
% Generate the network (represented by adjacency matrix and size vector) based on parameter set output from Latin Hypercube Design
[banksizedist,adjacency_matrix,probability_matrix] = netgen2(n_banks,LHD,numsamples);

% Compute global and local (i.e. centrality) measures) for each simulated network
[Networks,globalnetworkmeasures,localnetworkmeasures,Global_NM_table] = networkmeasures(adjacency_matrix,n_banks,numsamples,fig_output);
toc

save('simnet_250banks_1000networks.mat'); % Network data prior to cascading defaults model (for testing)

%--------------------------------------------------------------------------
%% Shock calibrations from paper
%%% Random shock (RS) : {1,5,10}/250 banks
%%% Targeted shock(TS): {1,5,10}/10 largest banks
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('simnet_250banks_1000networks.mat'); % Network data prior to cascading defaults model (for testing)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shockmagnitude  = 0.35;      % Fraction of external assets shocked
numshockedbanks = 5;       % Number of shocked banks       
shocktype       = 'targeted'; % Choose between 'random' or 'targeted' shock   
omega = 1.054;              % Calibrated semi-elasticity parameter in inverse demand function for external assets


if strcmp(shocktype,'random')
    targetgroup     = n_banks;
    shocklabel      = 'RS_';
elseif strcmp(shocktype,'targeted')
    targetgroup      = 0.04*n_banks; % = 10 banks
    shocklabel  = 'TS_';
end

shockpars       = [shockmagnitude numshockedbanks omega];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Toggle parameters related to high-level functioning of the algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dispoutput    = 0; % Boolean variable: Show output for propagation dynamics (slow: use for debugging when modifying code)
liquidityrisk = 0; % Boolean variable: Allow for indirect contagion via external asset firesales (set liquidityrisk == 1)

if liquidityrisk == 1
    shockpars(3) = omega;
    LRlabel      = '_LRon';
    LRval        = 'on';
else
    LRlabel = '_LRoff';
    LRval   = 'off';
end

%--------------------------------------------------------------------------
%% Cascading defaults model
%--------------------------------------------------------------------------

% Running the cascading defaults model on each simulated network
tic
for k = 1:numsamples
    ActiveBanks = 1:n_banks;
    FailedBanks = []; 
     
    banks   = struct;
    
    for i =1:n_banks
        banks(i).assets.total = banksizedist(i,k);
    end
    
    disp('=====================================================');
    fprintf(1,'%s shock of %d/%d banks, liquidity effects: %s\n',shocktype,numshockedbanks,targetgroup,LRval);
    disp('=====================================================');

    fprintf(1,'Network iteration: %d/%d\n',k,numsamples)
    
    BS_input = [LHD(k,7), LHD(k,8)];
    
    % NESTED FUNCTION: This constitutes the cascading defaults model that generates our results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Outer level: cascadingdefaults initiates bank balance sheets, sets up framework for averaging results across shocks.
    % Inner level: cascadealgorithm applies the shock and runs the contagion mechanism by which it propagates through the system
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [banks,cascademodelresults(k,:),BS_exp_vars(k,:),banks_CDM_results_av(:,:,k),banks_BS_results(:,:,k),...
        shockedbanks_initBS_av(k,:),shockedbanks_centrality_av(:,:,k),initshock_magnitude_av(k)] =...
        cascadingdefaults(banks,n_banks,ActiveBanks,FailedBanks,adjacency_matrix(:,:,k),probability_matrix(:,:,k),...
        localnetworkmeasures(:,:,k),BS_input,shockpars,shocktype,targetgroup,liquidityrisk,dispoutput);
    clc
end
toc
    
% Creating global (network-level) regression tables
CM_table = array2table(cascademodelresults);

CM_table.Properties.VariableNames = {'NumFailedBanks','Norm_FailedBanks','TotCapitalLoss','Norm_CapitalLoss','TotDepositLoss',...
    'Norm_DepositLoss','SimulationTime','d_assetprice_abs','d_assetprice_pct'};

BS_EV_table = array2table(BS_exp_vars);
BS_EV_table.Properties.VariableNames  = {'Init_Assets','Init_Capital','Norm_Init_Capital','Init_IB_exp','Norm_Init_IB_exp'};

% Creating local (bank level) regression tables
bank_results = cat(2,banks_CDM_results_av,banks_BS_results,localnetworkmeasures);

%--------------------------------------------------------------------------
%% Output results to .mat and .csv formats
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global analysis - organising output from cascadingdefaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (I) Centrality of shocked banks
for i = 1:numsamples
    SB_centrality(i,:) = shockedbanks_centrality_av(:,:,i);    
end

SBC_table = array2table(SB_centrality);
SBC_table.Properties.VariableNames = {'Indeg_cent_SB','Outdeg_cent_SB',...
                                      'Inclose_cent_SB','Outclose_cent_SB'...
                                      'Btwn_cent_SB','Pagerank_cent_SB'};

% (II) Initial shock magnitude                                 
initshock_magnitude_av = initshock_magnitude_av';
ISM_table = array2table(initshock_magnitude_av);
ISM_table.Properties.VariableNames = {'Shock_size'};

shockedbanks_initBS = shockedbanks_initBS_av;
IBS_table = array2table(shockedbanks_initBS);
IBS_table.Properties.VariableNames= {'Init_capital_SB','Init_assets_SB','Init_LR_SB','Init_ext_assets_SB','Init_IBB_SB','Init_IBL_SB'};

timestamp = datestr(datetime('today'));

Data       = [cascademodelresults BS_exp_vars globalnetworkmeasures,shockedbanks_initBS,SB_centrality,initshock_magnitude_av];
Data_table = [CM_table BS_EV_table  Global_NM_table IBS_table SBC_table ISM_table];

% cleaning workspace - 1
clearvars a_min a_ratio ActiveBank adjacency_matrix alpha beta banksizedist BSpars BS_input d gamma gamma_a i k LH_input n_banks...
    networkpars probability_matrix shockmagnitude shockpars shockpercentage shockspecs shocktype targetcriterion targetgroup...
    theta ActiveBanks FailedBanks shockedbanks_centrality_av

% (III) Saving output. Produces a string featuring: date produced, type of shock, number of shocked banks, presence of liquidity effects
save(strcat('output_',timestamp,'_',shocklabel,num2str(numshockedbanks),LRlabel,'.mat'));                            % Save workspace as .mat file
writetable(Data_table,strcat(data_output,'data_',timestamp,'_',shocklabel,num2str(numshockedbanks),LRlabel,'.csv')); % Output data to .csv file

% (IV) Extract 100 values from LHD for inclusion in paper
rand_LH_draws_ind = randsample(1:numsamples,100);

rand_LH_draws = LHD(rand_LH_draws_ind,:);
rand_LH_draws_tab = array2table(rand_LH_draws);
rand_LH_draws_tab.Properties.VariableNames = {'a_min','a_ratio','gamma_a','alpha','beta','d','theta','gamma'};

writetable(rand_LH_draws_tab,strcat(tab_output,'LH_draws.xlsx'));

% cleaning workspace - 2
clearvars timestamp shocklabel numshockedbanks rand_LH_draws numsamples
end
