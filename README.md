# Source Codes
The archive contains two .R files- SPNNC.R and auxillary.R . The SPNNC.R contains the source code of SPNNC. The main function to call is SPNNC(X,k1,k2,alpha), where X is a data frame containing the data points row-wise. The default value of the parameters are set as discussed in the paper. The function returns the adjacency matrix(G) of the nearest neighbours graph formed. The connected components of the graph corresponding to G give the different clusters.
\
The auxillary.R file can also be used here. Calling the function master.SPNNC(X,G) (where G is the adjacency matrix returned by SPNNC() of SPNNC.r when applied on X) returns a vector of labels, starting with 1, of the clusters as can be derived from G, in the order of the rows in X.

# Datasets
All the synthetic as well as real life datasets used in the experiments are available as .csv files. The last column of each file contains the ground truth.
