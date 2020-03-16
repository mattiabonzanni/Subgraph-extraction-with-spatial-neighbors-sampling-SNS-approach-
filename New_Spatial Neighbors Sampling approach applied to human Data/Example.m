%% SNS approach to any Matrix 
            %Analysis of a matrix
% To calculate SWP, CC and PL (real, lattice and random) of a given matrix.
    % INPUTS
% MatrixG= unidirected matrix to test;
    % OUTPUT:
% 1) The following outputs are structured in the FinalMatrixG table:
% Degree=degree of the graph G;
% SWPG=value of SWP of the Global graph 
% RegularCCG= Average Clustering Coef. value of the Lattice model of the Global graph
% NetCCG= Average Clustering Coef. value of the Global graph 
% RandCCG= Average Clustering Coef. value of the Random model of the global graph 
% RegularPLG= Average Path length value of the Lattice model of the Global graph  
% NetPLG= Average Path length value of the Global graph   
% RandPLG= Average Path length value of the Random model of the Global graph
            %To apply the SNS method to any input matrix
% Fix a spatial parameter rho and a seed node s  of index i(InitalSeedNode);
% Extract all the nodes with x coordinates in the range [xs- rho; xs+ rho] and y coordinates in the range [ys- rho; ys+ rho] and z coordinatez in the range [zs-rho;zs+rho]. An edge ij is preserved if and only if both nodei and nodej are retained in the subgraph; 
% Calculate SWPS, CCS and PLS of the subgraph created in the previous step;
% Repeat the process from the last two steps with the node of index i+1 until the index of the seed node is equal to LastSeedNode.
    %INPUTS
% MatrixG= unidirected matrix to test;
% rho= the spatial parameter used to select the subgraph area/volume;
% x= a vector composed by the x-coordinates assigned to each node (based on real world measuraments). The real coordinates were normalized in the range [0;1];
% y= a vector composed by the y-coordinates assigned to each node (based on real world measuraments).The real coordinates were normalized in the range [0;1];
% z= a vector composed by the z-coordinates assigned to each node (based on real world measuraments).The real coordinates were normalized in the range [0;1];
% InitialSeedNode= the first node used as seed node for the extraction. If not specified, 1 is assumed as default;
% LastSeedNode= the last node used as seed node for the extraction. If not specified, NG is assumed as default;
    %OUTPUT
% 2) The following values are structured in the FinalSubGraphs table:
% NodeID=the index of the node used as seed node;
% DegreeSubGraph= degree of the subgraph;
% FractionNodeSubGraph= the fraction of nodes exctract in the subgraph;
% SWPS= values of SWP of the subgraphs; 
% RegularCCS= Average Clustering Coef. values of the Lattice model of the subgraphs; 
% NetCCS= Average Clustering Coef. values of the subgraphs; 
% RandCCS= Average Clustering Coef. values of the Random model of the subgraphs; 
% RegularPLS= Average Path length values of the Lattice model of the subgraphs;  
% NetPLS= Average Path length values of the subgraphs;   
% RandPLS= Average Path length values of the Random model of the subgraphs;
    % Required codes:
% 1) avg_clus_matrix (written by Eric Bridgeford);
% 2) avg_path_matrix (written by Eric Bridgeford);
% 3) clustering_coef_matrix (code originally written by Mika Rubinov,UNSW, 2007-2010 and modified/written by Eric Bridgeford);
% 4) latmio_und_connected (Mika Rubinov, UNSW; Jonathan Power, WUSTL and Olaf Sporns, IU; Brain Connectivity Toolbox);
% 5) randmio_und_connected (Mika Rubinov, UNSW; Jonathan Power, WUSTL and Olaf Sporns, IU; Brain Connectivity Toolbox);
    % Written by Mattia Bonzanni

clear
clc
load Coactivation_matrix
MatrixG=Coactivation_matrix;
[FinalMatrix]=matrixAnalysis(MatrixG);

load x
load y
load z
rho=0.3;
[FinalSubGraph]=SNSmatrix(MatrixG, rho, x, y, z, 10, 12);