function [a,adjacency_matrix,probability_matrix] = netgen2(n_banks,X,numsamples)

%--------------------------------------------------------------------------

%%% Copyright (C) Nicolas K. Scholtes, 2017
%%% Distributed under GPL v3.0

%--------------------------------------------------------------------------

%-------------------------------------------------------------------------
%% Drawing node fitness distribution from a truncated power law
%-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation and preallocation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


a_min   = X(:,1);
a_max   = (a_min.*X(:,2));
gamma_a = (X(:,3));

alpha = X(:,4);
beta  = X(:,5);
d     = X(:,6);

% Monte Carlo simulation with $n_draws$ draws from the distribution using
% the inverse transform sampling method

n_draws = 100000;
n_notconnected = zeros(1,numsamples);
rnd_LB = 1;
rnd_UB = n_banks/5;

a_minvec = repmat(a_min,1,n_draws);
a_maxvec = repmat(a_max,1,n_draws);

% User-specified minimum degree of all nodes in the network

mindeg = 1;

a = zeros(n_banks,numsamples);

probability_matrix = zeros(n_banks);
adjacency_matrix = zeros(n_banks,n_banks,numsamples);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse transform sampling method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Drawing the bank size distribution');

a_diff = a_maxvec.^(1-gamma_a) - a_minvec.^(1-gamma_a);

u = rand(1,n_draws);

all_draws = (a_minvec.^(1-gamma_a)+ u.*(a_diff)).^(1./(1-gamma_a));

%all_draws = a_minvec.*(1-u).^(1./(1-gamma_a));
all_draws = all_draws';

% Logarithmic binning to obtain final fitness vector

for k = 1:numsamples
    [a(:,k),~] = lnbin(all_draws(:,k),n_banks);
end

a = sort(a,'descend');
a_max = max(a,[],1);

%-------------------------------------------------------------------------
%% Populating the adjacency matrix
%-------------------------------------------------------------------------

disp('Populating the adjacency matrix');

% Computing probability that two nodes form an undirected edge based on
% their relative fitness

rand_mat = rand(n_banks);

for k = 1:numsamples
    probability_matrix(:,:,k) = d(k)*((a(:,k)./a_max(k)).^alpha(k))*((a(:,k)./a_max(k)).^beta(k))';
    
    temp_PM = probability_matrix(:,:,k);
    temp_PM(1:n_banks+1:end) = 0;  
    probability_matrix(:,:,k) = temp_PM;
    
    connect = probability_matrix(:,:,k) >= rand_mat;
    adjacency_matrix(:,:,k) = double(connect);
  
end

for i=1:numsamples

    [S,~] = minimum_degree(adjacency_matrix(:,:,i),mindeg);
    
    if  ~S
        n_notconnected(i) = n_notconnected(i)+1;
    end
end

disp('Connecting isolated nodes');

for k=1:numsamples
        if ~minimum_degree(adjacency_matrix(:,:,k),mindeg)
           
            [~,isol_node_id_out] = ismember(adjacency_matrix(:,:,k),zeros(1,n_banks),'rows');
            [~,isol_node_id_in]  = ismember(adjacency_matrix(:,:,k)',zeros(1,n_banks),'rows');

            isol_node_out =  find(isol_node_id_out);
            isol_node_in  =  find(isol_node_id_in);
            
            isolated_nodes = intersect(isol_node_out,isol_node_in);
            
            n_isolatednodes = numel(isolated_nodes);
            connect_isolatednode = zeros(1,n_isolatednodes);
                        
            for i =1:n_isolatednodes
                connect_isolatednode(i) = round((rnd_UB-rnd_LB).*rand + rnd_LB);
                rewireprob = rand;
                if rewireprob>0.5
                    adjacency_matrix(isolated_nodes(i),connect_isolatednode(i),k)  = 1;
                else
                    adjacency_matrix(connect_isolatednode(i),isolated_nodes(i),k) = 1;
                end
            end
        end
    clearvars isolated_nodes isol_node_id_out isol_node_id_in isol_node_out isol_node_in 
end



end