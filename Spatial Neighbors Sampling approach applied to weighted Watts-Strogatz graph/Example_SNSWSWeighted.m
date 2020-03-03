%% SNS method applied to Weighted Watts-Strogatz graphs
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
% DegreeSubGraph=degree subgraph
% SWPG_values=values of SWP of the global graph 
% SWPS_values= values of SWP of the subgraphs 
% PercentageSubGraphsNodes= percentage of node to retain in each subGraph 
% RegularCCG= Average Clustering Coef. value of the Lattice model of the global graph
% NetCCG= Average Clustering Coef. value of the global graph 
% RandCCG= Average Clustering Coef. value of the Random model of the global graph 
% RegularCCS= Average Clustering Coef. value of the Lattice model of the subgraph  
% NetCCS= Average Clustering Coef. value of the subgraph   
% RandCCS= Average Clustering Coef. value of the Random model of the subgraph   
% RegularPLG= Average Path length value of the Lattice model of the global graph  
% NetPLG= Average Path length value of the global graph   
% RandPLG= Average Path length value of the Random model of the global graph   
% RegularPLS= Average Path length value of the Lattice model of the subgraph   
% NetPLS= Average Path length value of the subgraph 
% RandPLS= Average Path length value of the Random model of the subgraph
% The aformentioned outputs + the inputs (parameters) are structured in two tables:
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
%% Example
clc
NG=100;
k=15;
beta=0.05;          % rewiring probability
percentage=20;

[Final FinalMean]=SNSWSWeighted(NG, k, beta, percentage); 
