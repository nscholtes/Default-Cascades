function [reciprocity] = networkreciprocity(adj)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Code for computing reciprocity in unweighted, directed networks
% Based on Equation (2) in "Patterns of link reciprocity in directed networks"
% by Garlaschelli and Loffredo, 2008

% Copyright, Nicolas K. Scholtes, 2017
% Distributed under GPL v3.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Count number of bidirectional links for each node

L_bidir = 0;

for i = 1:size(adj,1)
    for j = 1:size(adj,1)
        if adj(i,j) == adj(j,i) && adj(i,j) == 1 && adj(j,i) == 1
            L_bidir = L_bidir+1;
        end
    end
end

% Compute traditional reciprocity measure given in Eq (1)

reciprocity_trad = L_bidir/sum(degrees_dir(adj));

% Compute simplified modified recipricity given in Eq (3)

reciprocity =  (reciprocity_trad - density_dir(adj))/(1- density_dir(adj));
            
end