function [S,num_no_incoming,num_no_outgoing] = minimum_degree(adj,mindeg)
 
% Check to see that all edges have a minimum user-specified degree
 
no_incoming = find(sum(adj) < mindeg);
no_outgoing = find(sum(adj,2) <mindeg);

num_no_incoming = numel(no_incoming);
num_no_outgoing = numel(no_outgoing);

% S is a binary parameter that is true only if there are NO isolated nodes and false if there exists at least one
 
if num_no_incoming ~= 0  && num_no_outgoing~=0 
    S= false;   
else
    S = true;
end
 
end