function [Networks,globalnetworkmeasures,localnetworkmeasures,Global_NM_table] = networkmeasures(adjacency_matrix,n_banks,numsamples,fig_output)

makeplots    = 0; % Select whether to output plots of various network measures
fig_output_N = strcat(fig_output,'Results-networks/');

node_nums = 1:n_banks;
node_ids  = cellfun(@num2str, num2cell(node_nums), 'UniformOutput', false);

mean_str = '$\bar{x}$ = '; var_str = '$s^{2}$ = ';

for k = 1:numsamples
    Network{k} = strcat('Network_',num2str(k));
end

bins = 50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preallocation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Global measures

density     = zeros(numsamples,1);
reciprocity = zeros(numsamples,1);

assortativity_oi = zeros(numsamples,1);
assortativity_io = zeros(numsamples,1);
assortativity_oo = zeros(numsamples,1);
assortativity_ii = zeros(numsamples,1);

distance_matrix = zeros(n_banks,n_banks,numsamples);
diameter        = zeros(numsamples,1);
avpathlength    = zeros(numsamples,1);
clustering_vec  = zeros(numsamples,n_banks);
avclustering    = zeros(numsamples,1);

out_degdist      = zeros(n_banks,numsamples);
in_degdist       = zeros(n_banks,numsamples);

indegree_centrality  = zeros(n_banks,numsamples);
outdegree_centrality = zeros(n_banks,numsamples);

incloseness_centrality  = zeros(n_banks,numsamples);
outcloseness_centrality = zeros(n_banks,numsamples);

betweenness_centrality = zeros(n_banks,numsamples);
pagerank_centrality    = zeros(n_banks,numsamples);

%--------------------------------------------------------------------------
%% GLOBAL NETWORK MEASURES
%--------------------------------------------------------------------------

% Density
for k=1:numsamples
    density(k) = density_dir(adjacency_matrix(:,:,k));
end

% reciprocity
for k = 1:numsamples
    reciprocity(k) = networkreciprocity(adjacency_matrix(:,:,k));
end

% Assortativity
for k=1:numsamples   
    assortativity_oi(k) = assortativity_bin(adjacency_matrix(:,:,k),1);
    assortativity_io(k) = assortativity_bin(adjacency_matrix(:,:,k),2);
    assortativity_oo(k) = assortativity_bin(adjacency_matrix(:,:,k),3);
    assortativity_ii(k) = assortativity_bin(adjacency_matrix(:,:,k),4);
end

assortativity = [assortativity_oi assortativity_io assortativity_oo assortativity_ii];

% Average path length
for k=1:numsamples    
    distance_matrix(:,:,k) = distance_bin(adjacency_matrix(:,:,k));
    [avpathlength(k),~,~,~,diameter(k)] = charpath(distance_matrix(:,:,k),0,0);
end

% Clustering coefficient

for k = 1:numsamples
    clustering_vec(k,:) = clustering_coef_bd(adjacency_matrix(:,:,k));
    avclustering(k)    = mean(clustering_vec(k,:));
end

% Degree distribution

for k = 1:numsamples
    Networks.(Network{k}) = digraph(adjacency_matrix(:,:,k),node_ids);
    
    out_degdist(:,k) = outdegree(Networks.(Network{k}));
    in_degdist(:,k)  = indegree(Networks.(Network{k}));
    
end

av_out_deg  = mean(out_degdist,2);   av_in_deg  = mean(in_degdist,2);
min_out_deg = min(out_degdist,[],2); max_in_deg = min(in_degdist,[],2);  
max_out_deg = max(out_degdist,[],2); max_in_deg = max(in_degdist,[],2); 

% Average degree

for k = 1:numsamples
    [indeg(k,:), outdeg(k,:),~] = degrees_dir(adjacency_matrix(:,:,k));
    
    av_indeg(k)  = mean(indeg(k,:));
    av_outdeg(k) = mean(outdeg(k,:));    
end

sys_av_indeg  = mean(av_indeg)
sys_av_outdeg = mean(av_outdeg);

sys_std_indeg = std(av_indeg)
sys_min_indeg = min(av_indeg)
sys_max_indeg = max(av_indeg)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output network measures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

globalnetworkmeasures = [density,diameter,reciprocity, assortativity,avpathlength,avclustering];

Global_NM_table = array2table(globalnetworkmeasures);

Global_NM_table.Properties.VariableNames = {'Density','Diameter','Reciprocity','OutIn_Assortativity','InOut_Assortativity','OutOut_Assortativity','InIn_Assortativity',...
    'Average_Path_Length','Average_Clustering'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting network measures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if makeplots == 1

% Density, average path length and clustering

mean_d = mean(density); mean_d_str = num2str(mean_d);
var_d  = var(density);  var_d_str  = num2str(var_d);

mean_apl = mean(avpathlength); mean_apl_str = num2str(mean_apl);
var_apl  = var(avpathlength);  var_apl_str  = num2str(var_apl);

mean_diam = mean(diameter); mean_diam_str = num2str(mean_diam);
var_diam  = var(diameter);  var_diam_str  = num2str(var_diam);

mean_cl = mean(avclustering); mean_cl_str = num2str(mean_cl);
var_cl  = var(avclustering);  var_cl_str  = num2str(var_cl);

figure
subplot(2,2,1)
    histogram(density)
    title('Network density','Interpreter','latex')
    dim = [0.3 0.55 0.3 0.35];
    str_d = {strcat(mean_str,mean_d_str),strcat(var_str,var_d_str)};
    annotation('textbox',dim,'String',str_d,'FitBoxToText','on','Interpreter','latex');
subplot(2,2,2)
    histogram(avpathlength)
    title('Average path length','Interpreter','latex')
    dim = [0.75 0.55 0.3 0.35];
    str_apl = {strcat(mean_str,mean_apl_str),strcat(var_str,var_apl_str)};
    annotation('textbox',dim,'String',str_apl,'FitBoxToText','on','Interpreter','latex'); 
subplot(2,2,3)
    histogram(diameter)
    title('Network diameter','Interpreter','latex')
    dim = [0.3 0.075 0.3 0.35];
    str_diam = {strcat(mean_str,mean_diam_str),strcat(var_str,var_diam_str)};
    annotation('textbox',dim,'String',str_diam,'FitBoxToText','on','Interpreter','latex'); 
subplot(2,2,4)
    histogram(avclustering)
    title('Average clustering','Interpreter','latex')
    dim = [0.75 0.075 0.3 0.35];
    str_cl = {strcat(mean_str,mean_cl_str),strcat(var_str,var_cl_str)};
    annotation('textbox',dim,'String',str_cl,'FitBoxToText','on','Interpreter','latex'); 
set(gcf,'renderer','painters');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-dpdf',strcat(fig_output_N,'Principalmeasures.pdf'));

% Reciprocity

mean_rec = mean(avclustering); mean_rec_str = num2str(mean_rec);
var_rec  = var(avclustering);  var_rec_str  = num2str(var_rec);

figure
histogram(reciprocity)
title('Reciprocity','Interpreter','latex')
dim = [0.3 0.55 0.3 0.35];
str_rec = {strcat(mean_str,mean_rec_str),strcat(var_str,var_rec_str)};
annotation('textbox',dim,'String',str_rec,'FitBoxToText','on','Interpreter','latex');
set(gcf,'renderer','painters');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-dpdf',strcat(fig_output_N,'Reciprocity.pdf'));


% Assortativity

mean_oi = mean(assortativity_oi); var_oi  = var(assortativity_oi); mean_oi_str = num2str(mean_oi); var_oi_str = num2str(var_oi);
mean_io = mean(assortativity_io); var_io  = var(assortativity_io); mean_io_str = num2str(mean_io); var_io_str = num2str(var_io);
mean_oo = mean(assortativity_oo); var_oo  = var(assortativity_oo); mean_oo_str = num2str(mean_oo); var_oo_str = num2str(var_oo);
mean_ii = mean(assortativity_ii); var_ii  = var(assortativity_ii); mean_ii_str = num2str(mean_ii); var_ii_str = num2str(var_ii);

mean_str = '$\bar{x}$ = '; var_str = '$s^{2}$ = ';

figure
subplot(2,2,1)
    histogram(assortativity(:,4))
    title('In-In Assortativity','Interpreter','latex')
    dim = [0.3 0.55 0.3 0.35];
    str = {strcat(mean_str,mean_ii_str),strcat(var_str,var_ii_str)};
    annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex');   
subplot(2,2,2)
    histogram(assortativity(:,2))
    title('In-Out Assortativity','Interpreter','latex')
    dim = [0.75 0.55 0.3 0.35];
    str = {strcat(mean_str,mean_io_str),strcat(var_str,var_io_str)};
    annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex');  
subplot(2,2,3)
    histogram(assortativity(:,1))
    title('Out-In Assortativity','Interpreter','latex')
    dim = [0.3 0.075 0.3 0.35];
    str = {strcat(mean_str,mean_oi_str),strcat(var_str,var_oi_str)};
    annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex');  
subplot(2,2,4)
    histogram(assortativity(:,3))
    title('Out-Out Assortativity','Interpreter','latex')
    dim = [0.75 0.075 0.3 0.35];
    str = {strcat(mean_str,mean_oo_str),strcat(var_str,var_oo_str)};
    annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex');  
set(gcf,'renderer','painters');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-dpdf',strcat(fig_output_N,'Assortativity.pdf'));

% In-degree distribution

[mean_id_counts, mean_id_bincenters] = hist(av_in_deg,100);[mean_out_counts, mean_out_bincenters] = hist(av_out_deg,100);

figure
[slope, ~,~, R2, S_mean_id] = logfit(mean_id_bincenters, mean_id_counts./sum(mean_id_counts),'powerlaw');

% Out-degree distribution

end

%--------------------------------------------------------------------------
%% LOCAL NETWORK MEASURES
%--------------------------------------------------------------------------

for k = 1:numsamples
    indegree_centrality(:,k)  = centrality(Networks.(Network{k}),'indegree');
    outdegree_centrality(:,k) = centrality(Networks.(Network{k}),'outdegree');
    
    incloseness_centrality(:,k)  = centrality(Networks.(Network{k}),'incloseness');
    outcloseness_centrality(:,k) = centrality(Networks.(Network{k}),'outcloseness');
    
    betweenness_centrality(:,k) = centrality(Networks.(Network{k}),'betweenness');
    pagerank_centrality(:,k)    = centrality(Networks.(Network{k}),'pagerank'); 
end

localnetworkmeasures = [reshape(indegree_centrality,[n_banks,1,numsamples])    reshape(outdegree_centrality,[n_banks,1,numsamples]) ...
                        reshape(incloseness_centrality,[n_banks,1,numsamples]) reshape(outcloseness_centrality,[n_banks,1,numsamples])...
                        reshape(betweenness_centrality,[n_banks,1,numsamples]) reshape(pagerank_centrality,[n_banks,1,numsamples])];
                    
%Local_NM_table = array2table(localnetworkmeasures);
%Local_NM_table.Properties.VariableNames = {'Indeg_ctrly','Outdeg_ctrly','Inclose_ctrly','Outclose_ctrly','Betweenness_ctrly','Pagerank_ctrly'};
end


