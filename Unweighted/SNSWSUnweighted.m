function [Final FinalMean]=SNSWSUnweighted(NG, k, beta, percentage)                         
%% SNS method applied to Unweighted Watts-Strogatz graphs
% 1) Given NG, k and beta (number nodes, average degree and re-wire probability, respectively), generate a Watts-Strogatz global network and calculate SWPG, CCG and PLG. ? were intentionally selected to achieve small-world networks; 
% 2) Fix the percentage of nodes NS to extract (percentage); 
% 3) Choose a seed node s of index i; 
% 4) Extract all the nodes with an index ? [i; i+NS]. The edges are maintained if and only if both nodes are retained in the subgraph; 
% 5) Calculate SWPS, CCS and PLS in the subgraph of dimension NS created in step3; 
% 6) Repeat the process from step3 to step5 with the node of index i+1 until the index of the seed node is equal to i-1 (to extract all the NG subgraphs of dimension NS).
   % INPUT:
% NG= number of nodes in the WS Global Graph;
% k= average degree in the WS Graph:
% beta= re-wiring probability in the WS Graph:0.005/0.01/0.05 are used to generate Small world networks. 
% percentage= percentage of node to retain in each subGraph. 
    % OUTPUT:
%DegreeSubGraph=degree subgraph
%SWPG_values=values of SWP of the global graph 
%SWPS_values= values of SWP of the subgraphs 
%PercentageSubGraphsNodes= percentage of node to retain in each subGraph 
%RegularCCG= Average Clustering Coef. value of the Lattice model of the global graph
%NetCCG= Average Clustering Coef. value of the global graph 
%RandCCG= Average Clustering Coef. value of the Random model of the global graph 
%RegularCCS= Average Clustering Coef. value of the Lattice model of the subgraph  
%NetCCS= Average Clustering Coef. value of the subgraph   
%RandCCS= Average Clustering Coef. value of the Random model of the subgraph   
%RegularPLG= Average Path length value of the Lattice model of the global graph  
%NetPLG= Average Path length value of the global graph   
%RandPLG= Average Path length value of the Random model of the global graph   
%RegularPLS= Average Path length value of the Lattice model of the subgraph   
%NetPLS= Average Path length value of the subgraph 
%RandPLS= Average Path length value of the Random model of the subgraph
% The aformentioned outputs + the inputs are structured in two tables:
% 1) Final=single values at each iteration;
% 2) FinalMean= mean values.
    %Required Code(s):
% 1) WattsStrogatz.m
% 2) avg_clus_matrix (written by Eric Bridgeford);
% 3) avg_path_matrix (written by Eric Bridgeford);
% 4) clustering_coef_matrix (code originally written by Mika Rubinov,UNSW, 2007-2010 and modified/written by Eric Bridgeford);
% 5) symm_matrix (written by Eric Bridgeford);
% 6) regular_matrix_generator (written by Eric Bridgeford);
% 7) randomize_matrix (written by Sarah F. Muldoon).
    % written by Mattia Bonzanni

%% Create Tables
nodeFraction=ceil(NG*percentage/100);               %to calculate the node's fraction
Parameters=array2table([NG k beta nodeFraction]);
TableSWP1 = cell(NG, 2);
TableSWP2 = cell(NG, 2);
TableVertex = cell(NG, 2);
TableNodeId = cell(NG, 2);
TableMeanDegreeSubH =cell(NG, 2);
TablePercentageSubGraphNodes=cell(NG, 2);
TableRegCC1=cell(NG, 2);
TableNetCC1=cell(NG, 2);
TableRandCC1=cell(NG, 2);
TableRegCC2=cell(NG, 2);
TableNetCC2=cell(NG, 2);
TableRandCC2=cell(NG, 2);
TableRegPL1=cell(NG, 2);
TableNetPL1=cell(NG, 2);
TableRandPL1=cell(NG, 2);
TableRegPL2=cell(NG, 2);
TableNetPL2=cell(NG, 2);
TableRandPL2=cell(NG, 2);
Tablenodes=cell(NG, 2);
%% Generate a Watts-Strogatz graph and a subgraph with defined distance across all the nodes

h=WattsStrogatz(NG,k,beta);                  % Generate a Watts-Strogatz graph - G(h)
degreeG=degree(h);                           % degree of global graph
MatrixGraph=full(adjacency(h));                     % Adj. Matrix of G(h)
%% SWP WS Graph
if sum(sum(MatrixGraph)) > 0    
bin_matrix = 0;
if strcmp('o','bin') == 1
   bin_matrix = 1;
   MatrixGraph = MatrixGraph > 0;
end
%check to see if matrix is symmeteric
symcheck=abs(MatrixGraph-MatrixGraph');
if sum(sum(symcheck)) > 0
    % adjust the input matrix to symmeterize
    disp('Input matrix is not symmetric. Symmetrizing.')
    W1 = symm_matrix(MatrixGraph, bin_matrix);
else
    W1=MatrixGraph;
end
%calculate the number of nodes
n1 = length(W1);  
%compute the weighted density of the network
dens_net1 = sum(sum(W1))/(max(max(W1))*n1*(n1-1));
%compute the average degree of the unweighted network, to give
%the approximate radius
numb_connections1 = length(find(W1>0));
avg_deg_unw1 = numb_connections1/n1;
avg_rad_unw1 = avg_deg_unw1/2;
avg_rad_eff1 = ceil(avg_rad_unw1);
%compute the regular and random matrix for the network W
W_reg1 = regular_matrix_generator(W1, avg_rad_eff1);        
W_rand1 = randomize_matrix(W1);
%compute all path length calculations for the network
reg_path1 = avg_path_matrix(1./W_reg1);     
rand_path1 = avg_path_matrix(1./W_rand1);   
net_path1 = avg_path_matrix(1./W1);          
p1 = (net_path1 - rand_path1);
if p1 < 0
    p1 = 0;
end
diff_path1 =  p1/ (reg_path1 - rand_path1);
if net_path1 == Inf || rand_path1 == Inf || reg_path1 == Inf
    diff_path1 = 1;
end
if diff_path1 > 1
    diff_path1 = 1;
end
%compute all clustering calculations for the network
reg_clus1 = avg_clus_matrix(W_reg1,'O');   
rand_clus1 = avg_clus_matrix(W_rand1,'O');
net_clus1 = avg_clus_matrix(W1,'O');
B1 = (reg_clus1 - net_clus1);
if B1 < 0
    B1 = 0;
end
diff_clus1 = B1 / (reg_clus1 - rand_clus1);
if isnan(reg_clus1) || isnan(rand_clus1) || isnan(net_clus1)
    diff_clus1 = 1;
end
if diff_clus1 > 1
    diff_clus1 = 1;
end
%calculate small world value, the root sum of the squares of
%diff path and diff clus
SWP1 = 1 - (sqrt(diff_clus1^2 + diff_path1^2)/sqrt(2));
delta_C1=diff_clus1;
delta_L1=diff_path1;
end
%% Extract all the SubGraphs and compute SWP
for i=1:NG
% Extract all the nodes with an index: [i; i+nodeFraction]. The edges are maintained if and only if both nodes are retained in the subgraph; 
Vertex=i;
nodeIDs=[Vertex:Vertex+nodeFraction];                         %to create a list of consecutive edges from Vertex to Vertex+Subnodes
indexesToReplace = nodeIDs > NG;                              %if Vertex+Subnodes>N, the next neighborn node is Vertex+Subnodes-N
nodeIDs(indexesToReplace) = nodeIDs(indexesToReplace)-NG;
sizeSubGraph=size(nodeIDs);
NodesSubGraph=sizeSubGraph(1,2);
ArraySubgraph=reshape(nodeIDs',1,[]);
SubH=subgraph(h,ArraySubgraph);                               %To generate a subGraph using the array "nodeIDs"
MeanDegree_SubH=mean(degree(SubH));                            %To calculate mean degree of each subgraphs
MatrixSub=full(adjacency(SubH));
                                                                    %%  SWP SubGraphs
if sum(sum(MatrixSub)) > 0    
bin_matrix = 0;
if strcmp('o','bin') == 1
   bin_matrix = 1;
   MatrixSub = MatrixSub > 0;
end
%check to see if matrix is symmeteric
symcheck1=abs(MatrixSub-MatrixSub');
if sum(sum(symcheck1)) > 0
    % adjust the input matrix to symmeterize
    disp('Input matrix is not symmetric. Symmetrizing.')
    W2 = symm_matrix(MatrixSub, bin_matrix);
else
    W2=MatrixSub;
end
%calculate the number of nodes
n2 = length(W2);  
%compute the weighted density of the network
dens_net2 = sum(sum(W2))/(max(max(W2))*n2*(n2-1));
%compute the average degree of the unweighted network, to give
%the approximate radius
numb_connections2 = length(find(W2>0));
avg_deg_unw2 = numb_connections2/n2;
avg_rad_unw2 = avg_deg_unw2/2;
avg_rad_eff2 = ceil(avg_rad_unw2);
%compute the regular and random matrix for the network W
W_reg2 = regular_matrix_generator(W2, avg_rad_eff2);        
W_rand2 = randomize_matrix(1./W2);
%compute all path length calculations for the network
reg_path2 = avg_path_matrix(1./W_reg2);      
rand_path2 = avg_path_matrix(1./W_rand2);    
net_path2 = avg_path_matrix(1./W2);          
p2 = (net_path2 - rand_path2);
if p2 < 0
    p2 = 0;
end
diff_path2 =  p2/ (reg_path2 - rand_path2);
if net_path2 == Inf || rand_path2 == Inf || reg_path2 == Inf
    diff_path2 = 1;
end
if diff_path2 > 1
    diff_path2 = 1;
end
%compute all clustering calculations for the network
reg_clus2 = avg_clus_matrix(W_reg2,'O');
rand_clus2 = avg_clus_matrix(W_rand2,'O');
net_clus2 = avg_clus_matrix(W2,'O');
B2 = (reg_clus2 - net_clus2);
if B2 < 0
    B2 = 0;
end   
diff_clus2 = B2 / (reg_clus2 - rand_clus2);
if isnan(reg_clus2) || isnan(rand_clus2) || isnan(net_clus2)
    diff_clus2 = 1;
end
if diff_clus2 > 1
    diff_clus2 = 1;
end
%calculate small world value, the root sum of the squares of
%diff path and diff clus
SWP2 = 1 - (sqrt(diff_clus2^2 + diff_path2^2)/sqrt(2));
delta_C2=diff_clus2;
delta_L2=diff_path2;
end
%% Store data
    TableSWP1(i, :) = {i SWP1};
    TableSWP2(i, :) = {i SWP2};
    TableVertex(i, :) = {i Vertex};
    TableNodeId (i, :)={i nodeIDs};
    TableMeanDegreeSubH (i, :)={i MeanDegree_SubH};
    TablePercentageSubGraphNodes (i, :)= {i percentage};
    TableRegCC1 (i, :) = {i reg_clus1};
    TableNetCC1 (i, :) = {i net_clus1};
    TableRandCC1 (i,:) = {i rand_clus1};
    TableRegCC2 (i, :) = {i reg_clus2};
    TableNetCC2 (i, :) = {i net_clus2};
    TableRandCC2 (i,:) = {i rand_clus2};
    TableRegPL1 (i, :) = {i reg_path1};
    TableNetPL1 (i, :) = {i net_path1};
    TableRandPL1 (i,:) = {i rand_path1};
    TableRegPL2 (i, :) = {i reg_path2};
    TableNetPL2 (i, :) = {i net_path2};
    TableRandPL2 (i,:) = {i rand_path2};
    Tableradius1 (i,:) = {i avg_rad_eff1};
    Tableradius2 (i,:) = {i avg_rad_eff2};
    Tablenodes (i,:) = {i nodeIDs};
end
                                                                 %% Compiles Tabeles
Table1 = cell2mat(TableSWP1);
Table2 = cell2mat(TableSWP2);
Table3= cell2mat(TableMeanDegreeSubH);
Table4=cell2mat(TablePercentageSubGraphNodes);
Table5=cell2mat(TableRegCC1);
Table6=cell2mat(TableNetCC1);
Table7=cell2mat(TableRandCC1);
Table8=cell2mat(TableRegCC2);
Table9=cell2mat(TableNetCC2);
Table10=cell2mat(TableRandCC2);
Table11=cell2mat(TableRegPL1);
Table12=cell2mat(TableNetPL1);
Table13=cell2mat(TableRandPL1);
Table14=cell2mat(TableRegPL2);
Table15=cell2mat(TableNetPL2);
Table16=cell2mat(TableRandPL2);
Table17=cell2mat(Tableradius1);
Table18=cell2mat(Tableradius2);
NodeIDs_values=Table1(: , 1);
SWPG_values=Table1(: , 2);
SWPS_values=Table2(: , 2);
DegreeSubGraph=Table3(: , 2);
PercentageSubGraphsNodes=Table4(: , 2);
RegularCCG=Table5(: , 2);
NetCCG=Table6(: , 2);
RandCCG=Table7(: , 2);
RegularCCS=Table8(: , 2);
NetCCS=Table9(: , 2);
RandCCS=Table10(: , 2);
RegularPLG=Table11(: , 2);
NetPLG=Table12(: , 2);
RandPLG=Table13(: , 2);
RegularPLS=Table14(: , 2);
NetPLS=Table15(: , 2);
RandPLS=Table16(: , 2);
                                                                        %% Output
S = [DegreeSubGraph SWPG_values SWPS_values PercentageSubGraphsNodes RegularCCG NetCCG RandCCG RegularCCS NetCCS RandCCS RegularPLG NetPLG RandPLG RegularPLS NetPLS RandPLS];
Final=array2table(full(S));
Final.Properties.VariableNames = {'DegreeSubGraph' 'SWPG_values' 'SWPS_values' 'PercentageSubGraphsNodes' 'RegularCCG' 'NetCCG' 'RandCCG' 'RegularCCS' 'NetCCS' 'RandCCS' 'RegularPLG' 'NetPLG' 'RandPLG' 'RegularPLS' 'NetPLS' 'RandPLS'};
FinalArray=full(S);
FinalMean=nanmean(FinalArray,1);
FinalMean=array2table(FinalMean);
FinalMean=cat(2,FinalMean, Parameters);
FinalMean.Properties.VariableNames = {'DegreeSubGraph' 'SWPG_values' 'SWPS_values' 'PercentageSubGraphsNodes' 'RegularCCG' 'NetCCG' 'RandCCG' 'RegularCCS' 'NetCCS' 'RandCCS' 'RegularPLG' 'NetPLG' 'RandPLG' 'RegularPLS' 'NetPLS' 'RandPLS' 'NG' 'k' 'beta' 'nodeFraction'};
clearvars -except degreeG W_reg1 W_rand1 Final FinalMean MatrixGraph MatrixSub net_clus1 net_clus2 net_path1 net_path2 rand_clus1 rand_clus2 rand_path1 rand_path2 reg_clus1 reg_clus2 reg_path1 reg_path2 SWP1 SWP2
z = msgbox('SNS Extraction from Unweighted Watts-Strogatz Graph Completed');
end