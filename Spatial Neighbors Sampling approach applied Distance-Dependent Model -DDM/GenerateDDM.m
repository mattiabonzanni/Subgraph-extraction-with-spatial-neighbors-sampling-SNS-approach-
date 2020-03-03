function [G MatrixG xcoord ycoord FinalGlobalGraph]=GenerateDDM(NG,edgeDensity,beta)                                                            
%% Generate the distance-dependent model (DDM) with small-world properties
% Given NG nodes, arrange them randomly in a 2D space. The 2D space has x- and y-coordinates range between 0 and 1;
% Assuming an inverse relationship between physical distance and edge strength (Muldoon, 2016), assign edge weights wij according to the euclidian distance dij between all the pairs of nodes as follows: w_ij=D_max-d_ij where Dmax=max{dij}. Each node will have NG-1 connections (except for the ij pair with distance d_ij equal to Dmax);
% Eliminate connections with an edge weight below a weight threshold wt. wt is defined by the user to achive the desired edge density. This step led to the spontaneous emergence of small-world networks;  
% Randomly re-wire each edge with probability ? (edge weight is retained). The inclusion of a random re-wiring guarantees that the network is not solely constructed as a function of physical distance, yet without imposing any additional rule. Moreover, the introduction of the random re-wiring led to a similar profile previously found in Watts et al.,1999. 
% Calculate SWPG, CCG and PLG.
    % INPUTS
% NG= number of nodes in the distance-dependent model (DDM) Global Graph;
% edgeDensity= the percentage of NG nodes attached to a given node in the Global graph.  If not specified, 8 is used (She et al., 2016);
% beta= re-wiring probability in the WS Graph. If not specified, 0.05 is used to generate Small world networks. 
    % OUTPUT:
% 1) G= the final DDM Global Graph;
% 2) MatrixG= the adjacency matrix of graph G;
% 3) xcoord= a vector composed by the random x-coordinates assigned to each node;
% 4) ycoord= a vector composed by the random y-coordinates assigned to each node;
% 5) The following outputs are structured in FinalGlobalGraph Table:
% Degree=degree of the graph G;
% wt=calculated weight thresholdto achieve the desired edge density;
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
% 4) latmio_und_connected (Mika Rubinov, UNSW; Jonathan Power, WUSTL and Olaf Sporns, IU);
% 5) randmio_und_connected (Mika Rubinov, UNSW; Jonathan Power, WUSTL and Olaf Sporns, IU);
    % Written by Mattia Bonzanni and Kimberly M. Bockley
%% Create DDN graph
xcoord=rand(1,NG);                                              % x-axis coordinates between [0;1]
ycoord=rand(1,NG);                                              % y-axis coordinates between [0;1]
nodeIds=1:NG;
distMatrix=zeros(NG);                                           % create an empty matrix
for i=1:1:NG             
    for j=2:1:NG         
        xDist= (xcoord(i)-xcoord(j))^2;                         % Calcuate squared distance between x values 
        yDist= (ycoord(i)-ycoord(j))^2;                         % Calcuate squared distance between y values
        distMatrix(i, j)= sqrt(xDist+yDist);                    % Calculate total distance using the pythagorean theorem
        distMatrix(j, i)= sqrt(xDist+yDist); 
    end
end
distMatrix=distMatrix-(max(distMatrix,[],'all'));               % to have closer nodes with higher weight
distMatrix=distMatrix./min(distMatrix,[],'all');                % to range every between [0;1]. It follows that we can have the threshold's range [0;1].
wt=0.1;                                                         % initial value of wt; the value will be updated until reaching the desired edge density
Degree=NG*(NG-1);
if edgeDensity==[]                                              % if edgeDensity is not decleared, 8 is set as default value, as previously reported in She et al.,2016 'Evaluating the Small-World-Ness of a Sampled Network: Functional Connectivity of Entorhinal-Hippocampal Circuitry' 
    edgeDensity=8;
end
while Degree>(2*edgeDensity*NG)/100                             % to identify the value of wt necessary to achieve the desired edge density
    distMatrix(distMatrix<=wt)=0;                               % to threshold the matrix in order to eliminate edges from nodes ij given the distance(ij)<threshold. 
    distMatrix = distMatrix - diag(diag(distMatrix));           % to eliminate self-loops (0 values on the diagonal)
    A=graph(distMatrix);                                    
    [s,t]=findedge(A);                                          % to extract the s and t table of the nodes pair
    weights=A.Edges.Weight(findedge(A,s,t));                    % to extract the edgeWeights table of the nodes pair
    NodeDegree=degree(A);                                       
    Degree=mean(NodeDegree);
    wt=wt+0.001;
    fprintf('wt=%.3f\n',wt);
end
wt=wt-0.001;
%% Re-wire edges keeping the Edge weights
if beta==[]                                                     % if beta is not decleared, 0.05 is set as default value. 
    beta=0.05;
end
source=rand(size(t,1),1);                                       % to generate random numbers between [0;1]; we need to generate as many as are the edges, thus I need to have size of table t
source(source<=beta)=1;                                         % to generate logic 1 values for numbers with probability equal or smaller than beta
source(source<1)=0;                                             % to generate logic 0 values for numbers with probability larger than beta
nodesToRewire=find(source);                                     % to find the index positions of the logic 1 values
howManyToRewire=size(nodesToRewire,1);                          % to find how many edges will be rewired
for q=1:howManyToRewire
    g=graph(s,t);                                               % creating the graph in the loop allows me to consider at every iteration the new targets for s in order to avoid multiple edges between pairs of nodes
    AllTargets=(1:NG);                                          % to list all the nodes
    position=nodesToRewire(q,1);                                % to substitute an edge in t using the position based on nodesToRewire
    nodeS=s(position,1);
    oldTargets=transpose(neighbors(g,nodeS));                   % to find the node in table s
    AllTargets(oldTargets)=[];                                  % to exclude old targets and thus generate a new list of avaible targets for rewiring
    newTargetidx=randi(size(AllTargets,2),1);                   % to generate random numbers to choose the new targets for rewiring from the reduced AllTargets array
    t(position,1)=AllTargets(1,newTargetidx);                   % to rewire in table t without multiple edges using a random value from the reduced AllTargets array  
end
G = graph(s,t,weights);
MatrixG=full(adjacency(G,'weighted'));                          % re-wired matrix - Global Graph
%% CC, PL and SWP analysis 
NetCCG=avg_clus_matrix(MatrixG, 'O'); 
NetPLG=avg_path_matrix(1./MatrixG);
RGlob=MatrixG;
Iter=200;
RandomMatrixG= randmio_und_connected(RGlob, Iter);
RegularMatrixG= latmio_und_connected(RGlob,Iter);
RandCCG= avg_clus_matrix(RandomMatrixG, 'O');
RegularCCG= avg_clus_matrix(RegularMatrixG, 'O');
RandPLG=avg_path_matrix(1./RandomMatrixG);
RegularPLG=avg_path_matrix(1./RegularMatrixG);
z = (NetPLG - RandPLG);
if z < 0
   z = 0;
end
diff_path =  z/ (RegularPLG - RandPLG);
if NetPLG == Inf || RandPLG == Inf || RegularPLG == Inf
    diff_path = 1;
end
if diff_path > 1
    diff_path = 1;
end

B = (RegularCCG - NetCCG);
if B < 0
    B = 0;
end
diff_clus = B / (RegularCCG - RandCCG);
if isnan(RegularCCG) || isnan(RandCCG) || isnan(NetCCG)
    diff_clus = 1;
end
if diff_clus > 1
    diff_clus = 1;
end

SWPG = 1 - (sqrt(diff_clus^2 + diff_path^2)/sqrt(2));

S = [NG beta Degree wt SWPG RegularCCG NetCCG RandCCG RegularPLG NetPLG RandPLG];
FinalGlobalGraph=array2table(full(S));
FinalGlobalGraph.Properties.VariableNames = {'NG' 'Beta' 'DegreeGlobalGraph' 'wt' 'SWPS_values' 'RegularCCS' 'NetCCS' 'RandCCS' 'RegularPLS' 'NetPLS' 'RandPLS'};
m1 = msgbox('Distance-Dependent Model (DDM) generated.');
end