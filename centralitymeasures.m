function [] = centralitymeasures(adjacency_matrix,Networks,numsamples,fig_output)

fig_output_N = strcat(fig_output,'Results-networks/');

%--------------------------------------------------------------------------
% Computing centrality measures
%--------------------------------------------------------------------------

% Degree centrality

% Betweenness centrality

% 

%--------------------------------------------------------------------------
% Plotting centrality measures
%--------------------------------------------------------------------------

figure
cc = centrality(Networks.Network_1,'betweenness');
cc(cc==0) = 1;
p=plot(Networks.Network_1,'Layout','force','NodeCData',cc,'MarkerSize',10*cc/max(cc));
labelnode(p,node_ids(1:10),node_ids(1:10));   
colormap jet
colorbar
title('Betweenness Centrality Scores - Unweighted')
set(gcf,'renderer','painters');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-dpdf',strcat(fig_output_N,'Betweenness_centrality.pdf'));

end