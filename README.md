# Subgraph-extraction-with-spatial-neighbors-sampling-SNS-approach-
MATLAB codes to extract spatial neighbors from graphs and raw data. 
Reference Paper:
Bonzanni M., Bockley K., Kaplan D.L. On the effect of neuronal spatial subsampling in small-world networks, Eur J Neurosci. 2020;00:1–14.  https://doi.org/10.1111/ejn.14937

  WHY:
The analysis of real world networks is biased by the current ability to measure just a subsample of the entire network.  
It is thus relevant to understand if the information gained in the subsamples can be extended to the global network. 
Here we showed how average clustering coefficient, average path length and small-world propensity scale when a spatial sampling is applied to small-world networks. 
This extraction mimics the measurement of physical neighbors by means of electrical and optical techniques, both used to study neuronal networks.

  HOW:
The spatial neighbors sampling (SNS) approach aims to retain spatial neighboring nodes in the subsample. 
The definition of spatial neighbors is either based on: 
1) node indexing (nodes whose indices are consecutive integers in a given interval as seen in the Watts-Strogatz models);
2) Euclidian distance (nodes with coordinates within a 2D surface/3D volume as seen in the distance-dependent model and human data, respectively). 

The extracted graph is an induced subgraph since it contains all the edges connecting pairs of retained nodes. 
Differently from node, link or snowball sampling methods, we did not use any topological information of the global network during the subgraph extraction. 
This is aimed to mimic the lack of information of the global network while optically selecting a single field-of-view during the study of neuronal network. 

There are four folders:
1) Spatial Neighbors Sampling approach applied to unweighted Watts-Strogatz graph

Generation of binary Watts-Strogatz graphs and subgraph extraction.

2) Spatial Neighbors Sampling approach applied to weighted Watts-Strogatz graph

Generation of weighted Watts-Strogatz graphs and subgraph extraction.

3) Spatial Neighbors Sampling approach applied Distance-Dependent Model -DDM

Generation of distance-dependent graphs and subgraph extraction.

4) New_Spatial Neighbors Sampling approach applied to human Data

Analysis of the human coactivation matrix (Crossley et al., 2013) and subgraph extraction. 

Each folder contains:
1) A function to generate a global graph and extract the subgraphs (Watts-Strogatz, unweighted and weighted) or two distinct functions to create/analyze a graph and extract subgraphs (distance-dependent model and human data);
2) An example file;
3) All the required codes. 

There is also a raw.m structure variable which contains:
1. Raw data for the binary and weighted Watts-Strogatz models;
2. Graphs and coordinates for the ditance-dependent model and the human dataset.

For any information, do not hesitate to contact me at mattia.bonzanni@tufts.edu or mattia.bonzanni@hotmail.it
