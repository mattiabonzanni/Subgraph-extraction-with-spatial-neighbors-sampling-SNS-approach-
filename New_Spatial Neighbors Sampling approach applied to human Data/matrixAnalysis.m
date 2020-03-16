function [FinalMatrix]=matrixAnalysis(MatrixG);
%% Analysis of a matrix
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
    % Required codes:
% 1) avg_clus_matrix (written by Eric Bridgeford);
% 2) avg_path_matrix (written by Eric Bridgeford);
% 3) clustering_coef_matrix (code originally written by Mika Rubinov,UNSW, 2007-2010 and modified/written by Eric Bridgeford);
% 4) latmio_und_connected (Mika Rubinov, UNSW; Jonathan Power, WUSTL and Olaf Sporns, IU; Brain Connectivity Toolbox);
% 5) randmio_und_connected (Mika Rubinov, UNSW; Jonathan Power, WUSTL and Olaf Sporns, IU,Brain Connectivity Toolbox);
    % Written by Mattia Bonzanni 

 %% To calcluate CC, PL and SWP     
NG=size(MatrixG,1);                                                 % number of nodes in the global graph
Degree=mean(degree(graph(MatrixG)));                                % degree of the global graph
NetCCG=avg_clus_matrix(MatrixG, 'O');                               % average clustering coeficient of the global graph
NetPLG=avg_path_matrix(1./MatrixG);                                 % average path length of the global graph
R=MatrixG;
Iter=200;                                                           % arbitrarly set at 200
RandomMatrixG= randmio_und_connected(R,Iter);
fprintf('   Randomization of the Graph is completed.\n');
LatticeMatrixG= latmio_und_connected(R,Iter);
fprintf('      Latticization of the Graph is completed.\n');
RandCCG= avg_clus_matrix(RandomMatrixG, 'O');
LattCCG= avg_clus_matrix(LatticeMatrixG, 'O');
RandPLG=avg_path_matrix(1./RandomMatrixG);
LatticePLG=avg_path_matrix(1./LatticeMatrixG);
z = (NetPLG - RandPLG);
if z < 0
    z = 0;
end
diff_pathG =  z/ (LatticePLG - RandPLG);
if diff_pathG > 1
    diff_pathG = 1;
end
B = (LattCCG - NetCCG);
if B < 0
    B = 0;
end
diff_clusG = B / (LattCCG - RandCCG);
if diff_clusG > 1
    diff_clusG = 1;
end
SWPG = 1 - (sqrt(diff_clusG^2 + diff_pathG^2)/sqrt(2));
%% Output
Final=table(NG, Degree, SWPG,LattCCG, NetCCG, RandCCG, LatticePLG, NetPLG, RandPLG);
FinalMatrix=array2table(Final);
FinalMatrixG.Properties.VariableNames = {'NG' 'DegreeGraph' 'SWPG' 'RegularCCG' 'NetCCG' 'RandCCG' 'RegularPLG' 'NetPLG' 'RandPLG'};
m4 = msgbox('Network Analysis Completed');
end