function [] = topsim()

%--------------------------------------------------------------------------
% MATLAB code for 'Default cascades and systemic risk on different multilayer interbank network topologies'
% by Nicolas K. Scholtes (2017)

% Copyright (C) Nicolas K. Scholtes, 2017
% Distributed under GPL 3.0

addpath(genpath('/Users/nscholte/Desktop/Research/Ch.3 - Systemic risk/MATLAB Codes'));


% Output figures to directory with LaTeX draft for automatic updating
fig_output  = '/Users/nscholte/Desktop/Research/Ch.3 - Systemic risk/Drafts/Figures/';
data_output = '/Users/nscholte/Desktop/Research/Ch.3 - Systemic risk/R Codes/';
tab_output  = '/Users/nscholte/Desktop/Research/Ch.3 - Systemic risk/Drafts/Tables/';

%--------------------------------------------------------------------------
%% CALIBRATION OF SIMULATION PARAMETERS
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bank size/fitness distribution parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_banks     = 250;        % Number of banks

a_min   = [10 50];         
a_ratio = [10,50];       
gamma_a = [2 3];         

alpha   = [0.5 1];
beta    = [0.5 1];
d       = [0.75 0.9];

networkpars = [a_min; a_ratio; gamma_a; alpha; beta; d];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial balance sheet weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta = 0.8;   % external asset/asset ratio
gamma = 0.1;   % capital/liabilites ratio

BSpars = [theta gamma];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shock calibrations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shockmagnitude  = 0.5;      % Fraction of external assets shocked
numshockedbanks = 1;

shocktype       = 'random'; % Choose between random or targeted shock
targetcriterion = 'betweenness_centrality'; % Location of targeted shock based on various criteria
                                            % Bank size or various node centrality measures
shockspecs = {shocktype,targetcriterion};                                           
shockpars  = [shockmagnitude numshockedbanks];

%--------------------------------------------------------------------------
% Latin hypercuble sampling
%--------------------------------------------------------------------------
numsamples = 1000;
[LHD] = createLH(networkpars,numsamples);
%--------------------------------------------------------------------------
% Generating the network
%--------------------------------------------------------------------------

tic
[banksizedist,adjacency_matrix,probability_matrix] = netgen2(n_banks,LHD,numsamples);
toc

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
    
    [~,cascademodelresults(k,:)] = cascadingdefaults(banks,n_banks,ActiveBanks,FailedBanks,...
        adjacency_matrix(:,:,k),probability_matrix(:,:,k),BSpars,shockpars,shockspecs);
end
toc

CM_table = array2table(cascademodelresults);

CM_table.Properties.VariableNames = {'NumFailedBanks','TotCapitalLoss','SimulationTime'}
load('simnet+casc_250banks_1000networks.mat');

%--------------------------------------------------------------------------
% Global network measures
%--------------------------------------------------------------------------

[Networks,globalnetworkmeasures,Global_NM_table] = networkmeasures(adjacency_matrix,n_banks,numsamples,fig_output);

%--------------------------------------------------------------------------
% Global regressions
%--------------------------------------------------------------------------

% Compute total capital initially in the system


Data       = [cascademodelresults globalnetworkmeasures];
Data_table = [CM_table Global_NM_table];

writetable(Data_table,strcat(data_output,'data.csv'));

globalregressions(Data);


end
