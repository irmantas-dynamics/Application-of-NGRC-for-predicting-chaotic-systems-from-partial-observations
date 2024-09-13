[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13756781.svg)](https://doi.org/10.5281/zenodo.13756781)

Code for reproducing the figures from the paper titled "Application of next-generation reservoir computing for predicting chaotic systems from partial observations".
See: https://doi.org/10.1103/PhysRevE.109.064215


The main structure of all programs is as follows:
1. Function `Generate_signal()` prepares a signal that later will be used for the training and the prediction stages. This function calls to function `eqsn()`, where equations of the dynamical system are defined;
2. The training signal `X0` and signal `Xpred` that will be used for the estimation of the prediction quality  are used as input for functions `Cheb_prediction()` or `Monom_prediction()`;
3. The functions `Cheb_prediction()` or `Monom_prediction()`:
   - Rescale training signal `X0`;
   - Generates file '*Monom_Rn_k.m*' or '*Cheb_Rn_k.m*', where k is the number of embedded dimensions. Functions  `Monom_Rn_k()` and `Cheb_Rn_k()` from a given signal `X` combine feature vector;
   - Compose an array of the feature vectors `R` within the call of function `features()`;
   - Finally, it calculates the output matrix `W`, that later can be used for the prediction.
