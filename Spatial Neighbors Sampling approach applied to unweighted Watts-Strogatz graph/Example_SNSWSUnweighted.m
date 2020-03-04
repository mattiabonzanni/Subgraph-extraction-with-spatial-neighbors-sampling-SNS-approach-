%% SNS method applied to Unweighted Watts-Strogatz (WS) graphs
% 1) Given NG, k and beta (number nodes, average degree and re-wire probability, respectively), generate a Watts-Strogatz global network and calculate SWPG, CCG and PLG. beta is intentionally selected to achieve small-world networks; 
% 2) Fix the percentage of nodes NS to extract (percentage); 
% 3) Choose a seed node s of index i; 
% 4) Extract all the nodes with an index in the range [i; i+NS]. The edges are maintained if and only if both nodes are retained in the subgraph; 
% 5) Calculate SWPS, CCS and PLS in the subgraph of dimension NS created in the previous step; 
% 6) Repeat the process (last three steps) with the node of index i+1 until the index of the seed node is equal to i-1 (to extract all the NG subgraphs of dimension NS).
   % INPUT:
% NG= number of nodes in the WS Global Graph;
% k= average degree in the WS Graph:
% beta= re-wiring probability in the WS Graph. 0.005/0.01/0.05 are used to generate Small world networks. 
% percentage= percentage of node to retain in each subGraph;
% InitialSeedNode= the first node used as seed node for the extraction. If not specified, 1 is assumed as default;
% LastSeedNode= the last node used as seed node for the extraction. If not specified, NG is assumed as default;
    % OUTPUT:
%DegreeSubGraph=degree subgraph;
%SWPG_values=values of SWP of the global graphs; 
%SWPS_values= values of SWP of the subgraphss 
%PercentageSubGraphsNodes= percentage of node retained in each subGraph; 
%RegularCCG= Average Clustering Coef. values of the Lattice model of the global graphs;
%NetCCG= Average Clustering Coef. values of the global graphs; 
%RandCCG= Average Clustering Coef. values of the Random model of the global graphs; 
%RegularCCS= Average Clustering Coef. value of the Lattice model of the subgraphs;  
%NetCCS= Average Clustering Coef. values of the subgraphs;   
%RandCCS= Average Clustering Coef. values of the Random model of the subgraphs;   
%RegularPLG= Average Path length values of the Lattice model of the global graphs;  
%NetPLG= Average Path length values of the global graphs;   
%RandPLG= Average Path length values of the Random model of the global graphs;   
%RegularPLS= Average Path length values of the Lattice model of the subgraphs   
%NetPLS= Average Path length values of the subgraphs; 
%RandPLS= Average Path length values of the Random model of the subgraphs;
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
%% Example
clear
clc
NG=100;             % number of nodes in the global WS graph
k=15;               % mean node degree 2*K
beta=0.05;          % rewiring probability
percentage=15;      % percentage of nodes to extract 

[Final FinalMean]=SNSWSUnweighted(NG, k, beta, percentage);         %If not specified, InitialSeedNode=1 and LasteSeedNode=NG.
