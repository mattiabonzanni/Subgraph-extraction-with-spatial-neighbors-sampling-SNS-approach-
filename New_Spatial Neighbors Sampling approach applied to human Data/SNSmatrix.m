function [FinalSubGraph]=SNSmatrix(MatrixG, rho, x, y, z, InitialSeedNode, LastSeedNode)
%% To apply the SNS method to any input matrix
% Fix a spatial parameter rho and a seed node s  of index i(InitalSeedNode);
% Extract all the nodes with x coordinates in the range [xs- rho; xs+ rho] and y coordinates in the range [ys- rho; ys+ rho] and z coordinatez in the range [zs-rho;zs+rho]. An edge ij is preserved if and only if both nodei and nodej are retained in the subgraph; 
% Calculate SWPS, CCS and PLS of the subgraph created in the previous step;
% Repeat the process from the last two steps with the node of index i+1 until the index of the seed node is equal to LastSeedNode.
    %INPUTS
% MatrixG= unidirected matrix to test;
% rho= the spatial parameter used to select the subgraph area/volume;
% x= a vector composed by the x-coordinates assigned to each node (based on real world measuraments).The real coordinates were normalized in the range [0;1];
% y= a vector composed by the y-coordinates assigned to each node (based on real world measuraments).The real coordinates were normalized in the range [0;1]);
% z= a vector composed by the z-coordinates assigned to each node (based on real world measuraments).The real coordinates were normalized in the range [0;1];
% InitialSeedNode= the first node used as seed node for the extraction. If not specified, 1 is assumed as default;
% LastSeedNode= the last node used as seed node for the extraction. If not specified, NG is assumed as default;
    %OUTPUT
% 1) The following values are structured in the FinalSubGraphs table:
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
%% To create empty tables
NG=size(MatrixG,1);                                                         % the size of the global graph
G=graph(MatrixG);
if nargin<6                                                                 % if not decleared, initial node is 1 and last node is NG. Namely, it extracts all the subgraphs. 
   InitialSeedNode=1;
   LastSeedNode=NG;
end
NumberSeeds=LastSeedNode-InitialSeedNode+1;
S=zeros(NumberSeeds,10);                                                    % preallocation matrix
%% SubGraph extraction and analyis
for i=InitialSeedNode:LastSeedNode                                          % to select the fraction of nodes as seed node         
    fprintf('Node:%d\n',i)
    load x;                                                                 % x-coordinates input from the Global Graph 
    load y;                                                                 % y-coordinates input from the Global Graph 
    load z;                                                                 % z-coordinates input from the Global Graph
    f=x(1,i);                                                               % to identify the x-coordinate of the ith node
    x1=f-rho;                                                               % to identify the lower limit of the x-coordinate range based on rho
    x2=f+rho;                                                               % to identify the upper limit of the x-coordinate range based on rho
    x(x<x1)=0;                                                              % to eliminate all the x-coordinates < lower limit
    x(x>x2)=0;                                                              % to eliminate all the x-coordinates > upper limit
    nodex=find(x);                                                          % to identify the index of all the nodes of which the x-coordinate is between the lower and upper limit
    q=y(1,i);                                                               % to identify the y-coordinate of the ith node
    y1=q-rho;                                                               % to identify the lower limit of the y-coordinate range based on rho
    y2=q+rho;                                                               % to identify the upper limit of the y-coordinate range based on rho
    y(y<y1)=0;                                                              % to eliminate all the y-coordinates < lower limit
    y(y>y2)=0;                                                              % to eliminate all the y-coordinates > upper limit
    nodey=find(y);                                                          % to identify the index of all the nodes of which the y-coordinate is between the lower and upper limit
    v=z(1,i);                                                               % to identify the y-coordinate of the ith node
    z1=v-rho;                                                               % to identify the lower limit of the y-coordinate range based on rho
    z2=v+rho;                                                               % to identify the upper limit of the y-coordinate range based on rho
    z(z<y1)=0;                                                              % to eliminate all the y-coordinates < lower limit
    z(z>y2)=0;                                                              % to eliminate all the y-coordinates > upper limit
    nodez=find(z);                                                          % to identify the index of all the nodes of which the z-coordinate is between the lower and upper limit
    intersectxy=intersect(nodex, nodey);                                    % to isolate all the nodes with the x-coordinate between the x-lower and x-upper limit AND with the Y-coordinate between the Y-lower and Y-upper limit
    nodesToExtract=intersect(intersectxy, nodez);                           % to isolate all the nodes with the x-coordinate between the x-lower and x-upper limit AND with the Y-coordinate between the Y-lower and Y-upper limit AND with the Z-coordinate between the Z-lower and Z-upper limit
    FractionNodeSubGraph=(size(nodesToExtract,2))/NG;                       % to calculate the number of nodes in the subGraph of seed node ith
    subH=subgraph(G, nodesToExtract);                                       % to extract the subgraph
    Subgraph=full(adjacency(subH,'weighted'));
    MeanDegree_Sub=mean(degree(subH));
    %% To calcluate CC, PL and SWP 
    NetCCsub=avg_clus_matrix(Subgraph, 'O'); 
    NetLsub=avg_path_matrix(1./Subgraph);
    Rsub=Subgraph;
    Iter=200;                                                               % arbitrarly set at 200
    if length(nodesToExtract)>3
        RandomMatrixsub= randmio_und_connected(Rsub, Iter);
        fprintf('   Randomization of the subgraph (seed node %i) is completed.\n',i);
        LatticeMatrixsub= latmio_und_connected(Rsub,Iter);
        fprintf('      Latticization of the subgraph (seed node %i) is completed.\n',i);
    else
        RandomMatrixsub=Rsub;
        LatticeMatrixsub=Rsub;
    end
    RandCCsub= avg_clus_matrix(RandomMatrixsub, 'O');
    LattCCsub= avg_clus_matrix(LatticeMatrixsub, 'O');
    RandLsub=avg_path_matrix(1./RandomMatrixsub);
    LatticeLsub=avg_path_matrix(1./LatticeMatrixsub);
    d = (NetLsub - RandLsub);
    if d < 0
        d = 0;
    end
    diff_path2 =  d/ (LatticeLsub - RandLsub);
    if NetLsub == Inf || RandLsub == Inf || LatticeLsub == Inf
        diff_path2 = 1;
    end
    if diff_path2 > 1
        diff_path2 = 1;
    end

    B = (LattCCsub - NetCCsub);
    if B < 0
        B = 0;
    end   
    diff_clus2 = B / (LattCCsub - RandCCsub);
    if isnan(LattCCsub) || isnan(RandCCsub) || isnan(NetCCsub)
        diff_clus2 = 1;
    end
    if diff_clus2 > 1
        diff_clus2 = 1;
    end

    SWP2 = 1 - (sqrt(diff_clus2^2 + diff_path2^2)/sqrt(2));
    subData = [i MeanDegree_Sub FractionNodeSubGraph SWP2 LattCCsub NetCCsub RandCCsub LatticeLsub NetLsub RandLsub];
    S(i-InitialSeedNode+1,:)=subData;
end
%% Output
FinalSubGraph=array2table(full(S));
FinalSubGraph.Properties.VariableNames = {'NodeID' 'DegreeSubGraph' 'FractionNodeSubGraph' 'SWPS' 'RegularCCS' 'NetCCS' 'RandCCS' 'RegularPLS' 'NetPLS' 'RandPLS'};
m3 = msgbox('SubGraphs Extraction Completed');
end