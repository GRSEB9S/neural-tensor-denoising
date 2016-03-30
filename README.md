# Tensor Denoising
Denoise neural signals using tensor decompositions.

# Data format
Data should be formatted as an array `Y` that is size `(N,T,C,R)` where `N` is the number of neurons, `T` is the number of time points, `C` is the number of conditions, and `R` is the total number of trials. The fourth index should be padded with `NaNs` in cases where each neuron-condition pair did not see the same number of trials. Further, trial ordering conveys no meaning.

# Main Files
* `tensorDenoiseGridSearchCV.m` -- Finds the optimal tensor rank via grid search and leave-one-out cross-validation.
* `tensorDenoiseCluster.m` -- A wrapper file to pass parameters through tensorDenoiseGridSearchCV.m.
* `tensorDenoiseERR.m` -- Computes the relative error of tensor Yhat relative to Y.
* `tensorDenoiseDat2Plot.m` -- Takes the output files from `tensorDenoiseCluster.m` and produces the cell array `errPlot`, which can be used to plot results. 