
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Map discontinuity based diagnosis for neighbor embedding methods

The package is an implementation based on the paper:

Liu, Z., Ma, R., Zhong, Y. (2024) Assessing and improving reliability of neighbor embedding methods: a map-continuity perspective. arXiv:2410.16608.

The package is tested under R version 4.2.1 on Mac OS, and requires the R packages: `Rfast` (version 2.1.0), `parellel` (version 4.2.1), `RtsneWithP`.

## Installation

To install the package from the github repository, use:

``` r
if(!require(devtools)) install.packages("devtools") # If not already installed
devtools::install_github("zhexuandliu/RtsneWithP")
devtools::install_github("zhexuandliu/neMDBD")
```

For Mac users, it is possible to encounter error when installing the package `RtsneWithP` if `gfortran` is not already installed. Please see https://github.com/zhexuandliu/MapContinuity-NE-Reliability/blob/main/RtsneWithP/gfortran-install.md for detailed instructions on the installation of `gfortran`.

## Usage

``` r
### load package and calculate embedding

# Load package
library(neMDBD)
library(RtsneWithP)
library(ggplot2)
library(Rfast)
library(MGMM)

# generate Gaussian mixture data as an example
set.seed(1)
X = MGMM::rGMM(300, d = 2, k = 3, means = list(c(2, 0), c(-2, 0), c(0,-4)), covs = diag(2))
label = factor(rownames(X))

# run t-SNE
perplexity = 75
PCA_x = prcomp(X)
tsne_out = RtsneWithP::Rtsne(X, perplexity = perplexity, theta = 0, 
                 max_iter = 1000, Y_init = PCA_x$x[, 1:2]) # set theta = 0 to run exact tSNE

ggplot() +
  geom_point(data = data.frame(x = X[, 1], y = X[, 2], label = label),aes(x = x, y = y, color = label)) + 
  ggtitle('Original data with labels')
```

![](tools/unnamed-chunk-3-1.png)<!-- -->

``` r

ggplot() + 
  geom_point(data = data.frame(tSNE1 = tsne_out$Y[, 1], tSNE2 = tsne_out$Y[, 2], label = label),aes(x = tSNE1, y = tSNE2, color = label)) + 
  ggtitle('Embedding with labels')
```

![](tools/unnamed-chunk-3-2.png)<!-- -->

### Perturbation score

The perturbation score is available for t-SNE algorithm.

#### Inputs

- `X`: Matrix. The original data to be visualized.
- `Y`: Matrix. The embedding of the original data.
- `perplexity`: Integer. The perplexity parameter used in the t-SNE
  algorithm to obtain `Y`.
- `length`: Numeric. Length of perturbation. This parameter must be
  specified by the user.
- `approx`: Integer. `0` for computing the exact perturbation scores.
  `1` implies reusing the original PCA result for approximation. `2`
  accelerates the algorithm by computing the approximation of the
  similarity matrix.
- `ind`: Vector. The indices of points for which perturbation scores
  will be calculated. (default: `NULL` (calculate for all points))
- `no.cores`: Integer. Number of cores for parallel computing. (default:
  NULL (no parallel computing))
- `initial_dims`: Integer. The number of dimensions that should be
  retained in the initial PCA step in t-SNE. (default: `50` (as in the
  `Rtsne` implementation))
- `pca_center`: Logical. Should data be centered before pca is applied?
  (default: `TRUE` (as in the `Rtsne` implementation))
- `pca_scale`: Logical. Should data be scaled before pca is applied?
  (default: `FALSE` (as in the `Rtsne` implementation))

#### Output

- `pscore`: Vector. Perturbation scores.

#### Tutorial

``` r
### calculate perturbation score for all points
### Set approx = 0 for exact perturbation score, approx = 1 or 2 for acceleration
pscore = perturbation_score_compute(X, tsne_out$Y, perplexity, length = 0.5, approx = 0) # took 38.076s on a MacBook Air (M2 chip)
```

``` r
ggplot() +
  geom_point(data = data.frame(tSNE1 = tsne_out$Y[, 1], tSNE2 = tsne_out$Y[, 2], score = pscore),
             aes(x = tSNE1, y = tSNE2, color = score)) +
  viridis::scale_color_viridis(direction = 1, name = "Perturbation\nScore") + 
  ggtitle('Embedding with perturbation score')
```

![](tools/unnamed-chunk-5-1.png)<!-- -->

``` r
ggplot() +
  geom_point(data = data.frame(x = X[, 1], y = X[, 2], score = pscore),
             aes(x = x, y = y, color = score)) +
  viridis::scale_color_viridis(direction = 1, name = "Perturbation\nScore") + 
  ggtitle('Original data with perturbation score')
```

![](tools/unnamed-chunk-5-2.png)<!-- -->

##### Approximation method 1

``` r
# Set approx = 1 for acceleration by reusing the original PCA result
### The perturbation scores by the approximation method 1 are close to the exact perturbation scores.
pscore1 = perturbation_score_compute(X, tsne_out$Y, perplexity, length = 0.5, approx = 1) # took 36.665s on a MacBook Air (M2 chip) 
# The improvement is more significant for larger datasets as the PCA is relatively fast for small datasets.
ggplot() +
  geom_point(data = data.frame(approx0 = pscore, approx1 = pscore1),
             aes(x = approx0, y = approx1))
```

![](tools/unnamed-chunk-6-1.png)<!-- -->

##### Approximation method 2

``` r
### Set approx = 2 for acceleration by approximating the similarity matrix
### The perturbation scores by the approximation method 2 are close to the exact perturbation scores.
pscore2 = perturbation_score_compute(X, tsne_out$Y, perplexity, length = 0.5, approx = 2) # took 32.160s on a MacBook Air (M2 chip) 
# The improvement is more significant for larger datasets as the similarity matrix computation is relatively fast for small datasets.
ggplot() +
  geom_point(data = data.frame(approx0 = pscore, approx2 = pscore2),
             aes(x = approx0, y = approx2))
```

![](tools/unnamed-chunk-7-1.png)<!-- -->

##### Parallel computing

``` r
### set the number of cores `no.cores` for parallel computing
library(parallel)
pscore = perturbation_score_compute(X, tsne_out$Y, perplexity, length = 0.5, approx = 0, no.cores = 8) # took 13.355s on a MacBook Air (M2 chip)
```

##### Pre-screening points

``` r
### Pre-screening of points on the peripheries using 'dbscan'
### Most of the points with large perturbation scores are identified.
library(dbscan)
ind = which(!dbscan::is.corepoint(tsne_out$Y, eps = 1, minPts = 10))
pscore = rep(NA, dim(X)[1])
pscore[ind] = perturbation_score_compute(X, tsne_out$Y, perplexity, length = 0.5, approx = 2, ind = ind, no.cores = 8) # took 2.499s on a MacBook Air (M2 chip)
```

``` r
ggplot() +
  geom_point(data = data.frame(tSNE1 = tsne_out$Y[, 1], tSNE2 = tsne_out$Y[, 2], score = pscore),
             aes(x = tSNE1, y = tSNE2, color = score)) +
  viridis::scale_color_viridis(direction = 1, name = "Perturbation\nScore") + 
  ggtitle('Embedding with perturbation score (prescreened)')
```

![](tools/unnamed-chunk-10-1.png)<!-- -->

``` r
ggplot() +
  geom_point(data = data.frame(x = X[, 1], y = X[, 2], score = pscore),
             aes(x = x, y = y, color = score)) +
  viridis::scale_color_viridis(direction = 1, name = "Perturbation\nScore") + 
  ggtitle('Original data with perturbation score (prescreened)')
```

![](tools/unnamed-chunk-10-2.png)<!-- -->

### Singularity score

The singularity score is available for t-SNE algorithm.

#### Inputs

- `Y`: Matrix. The embedding of the original data.
- `P`: Matrix. The similarity matrix output by the t-SNE algorithm.

#### Output

- `sscore`: Vector. Singularity scores.

#### Tutorial

First, we should examine the gradient of the t-SNE loss function to
determine if a local minimum has been reached. Small gradient values
indicate convergence to a local minimum.

``` r
### check whether local minimum is reached
gradient_compute(tsne_out$Y, tsne_out$P)[1:10]
#>  [1] -9.002258e-17 -1.692569e-17 -1.887777e-17 -8.305646e-17  1.549520e-17
#>  [6]  1.774534e-18 -4.506946e-17  2.804973e-17  5.381179e-17 -5.280909e-17
```

After reaching the local minimum, we can calculate the singularity score
for fracture-inducing discontinuity diagnosis.

``` r
### calculate singularity score

# compute the singularity score
sscore = singularity_score_compute(tsne_out$Y, tsne_out$P) # took 0.055s on a MacBook Air (M2 chip)

# plot singularity score
ggplot() +
  geom_point(data = data.frame(tSNE1 = tsne_out$Y[, 1], tSNE2 = tsne_out$Y[, 2], sscore = sscore), aes(x = tSNE1, y = tSNE2, color = sscore)) +
  viridis::scale_color_viridis(direction = 1, trans = 'log10', name = "Singularity\nScore")
```

![](tools/unnamed-chunk-12-1.png)<!-- -->

Singularity scores can also guide the selection of perplexity. We select
the perplexity as the elbow point in the plot of the mean of top 5%
singularity scores versus perplexity.

``` r
### utilize singularity score to choose perplexity
### Select the perplexity as the elbow point in the plot of mean of top 5% singularity scores versus perplexity.

# calculate the singularity score for each perplexity candidate
perplexity_candidates = c(seq(5,90,5))
sscore_mat= matrix(NA, nrow = dim(X)[1], ncol = length(perplexity_candidates))
for (i in c(1:length(perplexity_candidates))){
  tsne_out_i = RtsneWithP::Rtsne(X, perplexity = perplexity_candidates[i], theta = 0, max_iter = 1000, Y_init = PCA_x$x[, 1:2])
  sscore_mat[,i] = singularity_score_compute(tsne_out_i$Y, tsne_out_i$P)
}

# plot the mean of top 5% singularity scores versus perplexity
# choose the elbow point as the perplexity to use
ggplot(data = data.frame(
  mean = apply(sscore_mat, 2, 
               function(s_score) {mean(s_score[s_score > quantile(s_score, 0.95)], na.rm = TRUE)}),
  perplexity = perplexity_candidates
)) +
  geom_point(aes(x = perplexity, y = mean))
```

![](tools/unnamed-chunk-13-1.png)<!-- -->

# Details

This R package implements functions to calculate singularity score and
perturbation score, to diagnose fracture-inducing and
overconfidence-inducing discontinuities in \[4\].

# References

\[1\] <https://github.com/jkrijthe/Rtsne>

\[2\] <https://github.com/jdonaldson/rtsne>

\[3\] L.J.P. van der Maaten and G.E. Hinton. “Visualizing
High-Dimensional Data Using t-SNE.” Journal of Machine Learning Research
9(Nov):2579-2605, 2008.

\[4\] Z. Liu, R. Ma, Y. Zhong. Assessing and improving reliability of
neighbor embedding methods: a map-continuity perspective.
arXiv:2410.16608, 2024.
