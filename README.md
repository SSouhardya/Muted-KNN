# Source Codes
The .R file SNNC.R contains the source code of SNNC. The main function to call is `SNNC(X,k,alpha,border_times,mode,cutoff)`, where X is a data frame containing the data points row-wise. The default value of the parameters are set as discussed in the paper.<br/><br/>
The default value of 'mode' is 1, which uses the Mahalanobis distance for outlier detection. To use the Distance Spike Method, set this value to 2. If 'mode' has value 2, then 'cutoff' needs to be specified, for the cutoff to be used in the DSS method. If the value is unspecified, the function creates a DSS plot that allows the user to choose a cutoff visually. In this case, the function doesn't resturn anything.<br/><br/>
With valid inputs, the function returns the cluster-labels (starting wth 1) for each point, as decided by the SNNC algorithm.


# Datasets
All the synthetic as well as real life datasets used in the experiments are available as .csv files. The last column of each file contains the ground truth.
