# Source Codes
The .R file SNNC.R contains the source code of SNNC. The main function to call is SNNC(X,k1,k2,alpha), where X is a data frame containing the data points row-wise. The default value of the parameters are set as discussed in the paper. The function returns the cluster-labels (startimg wth 1) for each point as decided by the SNNC algorithm. The program uses the package 'igraph' for searching connected components. The package should be installed. The first line of the code is a comment and contains the code for loading the same package in R.


# Datasets
All the synthetic as well as real life datasets used in the experiments are available as .csv files. The last column of each file contains the ground truth.
