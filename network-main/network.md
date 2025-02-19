Psychometric Network Analysis
================

## Setting Up R

- **R version:** R-4.4.2 for Windows (<https://cran.r-project.org/>)
- **Rtools version:** Rtools44 for Windows
  (<https://cran.r-project.org/bin/windows/Rtools/rtools44/rtools.html>)
- **RStudio version:** RStudio 2024.12.0+467
- **More information on EGAnet** <https://r-ega.net/>

To download the dataset:

<https://raw.githubusercontent.com/okanbulut/network/master/scr_canada.csv>

To download the R markdown file:

<https://raw.githubusercontent.com/okanbulut/network/master/network.Rmd>

To download the entire repository:

<https://github.com/okanbulut/network/archive/refs/heads/main.zip>

``` r
install.packages(c("dplyr", "ggcorrplot", "psych", "bootnet", "psychonetrics", "EGAnet", "networktree"))

library("dplyr")
library("ggcorrplot")
library("psych")
library("bootnet")
library("psychonetrics")
library("EGAnet")
library("networktree")
```

# Student Confident in Reading (SCR) in PIRLS 2021

The dataset comes from the PIRLS website
(<https://pirls2021.org/data/>).

``` r
# Students from four provinces in Canada
scr_canada <- read.csv("scr_canada.csv", header = TRUE)

# First 6 rows of the data
head(scr_canada)
```

    ##   IDCNTRY IDPOP IDGRADER IDGRADE IDSCHOOL IDCLASS   IDSTUD ITSEX ITADMINI
    ## 1    9134     1        2       4     5001  500101 50010101     2        2
    ## 2    9134     1        2       4     5001  500101 50010102     2        2
    ## 3    9134     1        2       4     5001  500101 50010103     1        2
    ## 4    9134     1        2       4     5001  500101 50010104     2        2
    ## 5    9134     1        2       4     5001  500101 50010105     2        2
    ## 6    9134     1        2       4     5001  500101 50010106     2        2
    ##   ITLANG_SA ITLANG_SQ IDBOOK ASBR08A ASBR08B ASBR08C ASBR08D ASBR08E ASBR08F
    ## 1         1         1      1       1       1       3       3       4       3
    ## 2         1         1     31       2       1       3       3       3       3
    ## 3         1         1      2       2       1       2       9       2       1
    ## 4         1         1     18       1       1       3       3       4       4
    ## 5         1         1     32       2       2       2       4       4       4
    ## 6         1         1     10       1       1       1       3       4       4
    ##     HOUWGT   TOTWGT   SENWGT JKREP JKZONE ASRREA01 ASRREA02 ASRREA03 ASRREA04
    ## 1 1.486722 22.28132 0.246146     1     70 599.3167 564.9619 577.9405 586.3821
    ## 2 1.486722 22.28132 0.246146     1     70 616.6323 589.4565 637.1342 625.7600
    ## 3 1.486722 22.28132 0.246146     1     70 548.8263 575.0935 565.4134 573.0199
    ## 4 1.486722 22.28132 0.246146     1     70 619.5068 669.4913 617.9325 607.4834
    ## 5 1.486722 22.28132 0.246146     1     70 586.4464 654.7330 621.8628 631.6019
    ## 6 1.486722 22.28132 0.246146     1     70 548.9577 620.0773 579.7630 561.6681
    ##   ASRREA05 idbid
    ## 1 590.1711   RR5
    ## 2 660.9687   RR5
    ## 3 578.2300   RR5
    ## 4 628.2404   RR5
    ## 5 600.9162   RR5
    ## 6 605.2696   RR5

``` r
# Variable names
names(scr_canada)
```

    ##  [1] "IDCNTRY"   "IDPOP"     "IDGRADER"  "IDGRADE"   "IDSCHOOL"  "IDCLASS"  
    ##  [7] "IDSTUD"    "ITSEX"     "ITADMINI"  "ITLANG_SA" "ITLANG_SQ" "IDBOOK"   
    ## [13] "ASBR08A"   "ASBR08B"   "ASBR08C"   "ASBR08D"   "ASBR08E"   "ASBR08F"  
    ## [19] "HOUWGT"    "TOTWGT"    "SENWGT"    "JKREP"     "JKZONE"    "ASRREA01" 
    ## [25] "ASRREA02"  "ASRREA03"  "ASRREA04"  "ASRREA05"  "idbid"

``` r
# Select only the SCR items
scr <- dplyr::select(scr_canada, starts_with("ASBR08"))

# Remove students with no valid responses and recode 9 as missing
scr %<>%
  dplyr::filter_all(all_vars(!is.na(.))) %>%
  dplyr::mutate_all(~na_if(., 9)) %>%
  as.data.frame()
```

Let’s check the correlation matrix of the items in SCR.

``` r
# Save the correlation matrix
cormat <- psych::polychoric(x = scr)$rho

# Correlation matrix plot
ggcorrplot::ggcorrplot(corr = cormat, # correlation matrix
                       type = "lower", # print only the lower part of the matrix
                       hc.order = TRUE, # hierarchical clustering
                       show.diag = TRUE, # show the diagonal values of 1
                       lab = TRUE, # add correlation values as labels
                       lab_size = 3) # Size of the labels
```

![](network_files/figure-gfm/ch3-1.png)<!-- -->

The correlation matrix above suggests that item wording effects might be
presented in the SCR scale. Let’s keep the original correlation matrix
and start analyzing the data using a Gaussian Graphical Model (GGM).

# Gaussian Graphical Model

``` r
network1 <- bootnet::estimateNetwork(
  data = scr, 
  corMethod = "cor_auto", # for polychoric and polyserial correlations
  default = "EBICglasso", # for estimating GGM with gLASSO and EBIC
  tuning = 0.5 # EBIC tuning parameter; set to zero for BIC model selection
)

# Print the estimated network
print(network1)
```

    ## 
    ## === Estimated network ===
    ## Number of nodes: 6 
    ## Number of non-zero edges: 15 / 15 
    ## Mean weight: 0.1088068 
    ## Network stored in network1$graph 
    ##  
    ## Default set used: EBICglasso 
    ##  
    ## Use plot(network1) to plot estimated network 
    ## Use bootnet(network1) to bootstrap edge weights and centrality indices 
    ## 
    ## Relevant references:
    ## 
    ##      Friedman, J. H., Hastie, T., & Tibshirani, R. (2008). Sparse inverse covariance estimation with the graphical lasso. Biostatistics, 9 (3), 432-441.
    ##  Foygel, R., & Drton, M. (2010). Extended Bayesian information criteria for Gaussian graphical models. 
    ##  Friedman, J. H., Hastie, T., & Tibshirani, R. (2014). glasso: Graphical lasso estimation of gaussian graphical models. Retrieved from https://CRAN.R-project.org/package=glasso
    ##  Epskamp, S., Cramer, A., Waldorp, L., Schmittmann, V. D., & Borsboom, D. (2012). qgraph: Network visualizations of relationships in psychometric data. Journal of Statistical Software, 48 (1), 1-18.
    ##  Epskamp, S., Borsboom, D., & Fried, E. I. (2016). Estimating psychological networks and their accuracy: a tutorial paper. arXiv preprint, arXiv:1604.08462.

``` r
# View the estimated network
plot(network1) 
```

![](network_files/figure-gfm/ch4-1.png)<!-- -->

Let’s also try the same model using psychonetrics.

``` r
# Save item names
obsvars <- colnames(scr)

network2 <- psychonetrics::ggm(scr, vars = obsvars) %>%
  psychonetrics::runmodel()

# Check out model fit
network2 %>% psychonetrics::fit()
```

    ##            Measure     Value
    ##               logl -84797.63
    ##  unrestricted.logl -84797.63
    ##      baseline.logl -97894.55
    ##               nvar         6
    ##               nobs        27
    ##               npar        27
    ##                 df       ~ 0
    ##          objective      2.59
    ##              chisq       ~ 0
    ##             pvalue         1
    ##     baseline.chisq  26193.85
    ##      baseline.npar        12
    ##        baseline.df        15
    ##    baseline.pvalue       ~ 0
    ##                nfi         1
    ##               pnfi       ~ 0
    ##                tli          
    ##               nnfi         1
    ##                rfi          
    ##                ifi         1
    ##                rni         1
    ##                cfi         1
    ##              rmsea          
    ##     rmsea.ci.lower       ~ 0
    ##     rmsea.ci.upper       ~ 0
    ##       rmsea.pvalue       ~ 0
    ##             aic.ll 169649.26
    ##            aic.ll2 169649.38
    ##              aic.x       ~ 0
    ##             aic.x2        54
    ##                bic 169849.87
    ##               bic2 169764.07
    ##            ebic.25 169898.25
    ##             ebic.5 169946.63
    ##            ebic.75 169985.33
    ##              ebic1 170043.38

We can prune this model to remove insignificant edges.

``` r
network3 <- psychonetrics::ggm(scr, vars = obsvars) %>%
  psychonetrics::runmodel() %>%
  psychonetrics::prune(adjust = "BH", alpha = 0.05)

# View the model parameters
network3 %>% psychonetrics::parameters()
```

    ## 
    ##  Parameters for group fullsample
    ##  -  mu  
    ##     var1 op var2  est     se        p row col par
    ##  ASBR08A ~1      1.49 0.0063 < 0.0001   1   1   1
    ##  ASBR08B ~1      1.53 0.0066 < 0.0001   2   1   2
    ##  ASBR08C ~1      2.58 0.0094 < 0.0001   3   1   3
    ##  ASBR08D ~1      3.14 0.0091 < 0.0001   4   1   4
    ##  ASBR08E ~1      3.31 0.0089 < 0.0001   5   1   5
    ##  ASBR08F ~1      3.40 0.0084 < 0.0001   6   1   6
    ## 
    ##  -  omega (symmetric) 
    ##     var1 op    var2    est     se        p row col par
    ##  ASBR08B -- ASBR08A   0.52 0.0065 < 0.0001   2   1   7
    ##  ASBR08C -- ASBR08A  0.042 0.0089 < 0.0001   3   1   8
    ##  ASBR08D -- ASBR08A -0.032 0.0090  0.00030   4   1   9
    ##  ASBR08E -- ASBR08A -0.037 0.0089 < 0.0001   5   1  10
    ##  ASBR08F -- ASBR08A  -0.15 0.0088 < 0.0001   6   1  11
    ##  ASBR08C -- ASBR08B -0.077 0.0089 < 0.0001   3   2  12
    ##  ASBR08D -- ASBR08B -0.082 0.0089 < 0.0001   4   2  13
    ##  ASBR08E -- ASBR08B -0.078 0.0089 < 0.0001   5   2  14
    ##  ASBR08F -- ASBR08B -0.076 0.0089 < 0.0001   6   2  15
    ##  ASBR08D -- ASBR08C   0.33 0.0080 < 0.0001   4   3  16
    ##  ASBR08E -- ASBR08C  0.076 0.0089 < 0.0001   5   3  17
    ##  ASBR08F -- ASBR08C  0.052 0.0089 < 0.0001   6   3  18
    ##  ASBR08E -- ASBR08D   0.37 0.0077 < 0.0001   5   4  19
    ##  ASBR08F -- ASBR08D   0.22 0.0085 < 0.0001   6   4  20
    ##  ASBR08F -- ASBR08E   0.34 0.0079 < 0.0001   6   5  21
    ## 
    ##  -  delta (diagonal) 
    ##     var1  op    var2  est     se        p row col par
    ##  ASBR08A ~/~ ASBR08A 0.54 0.0034 < 0.0001   1   1  22
    ##  ASBR08B ~/~ ASBR08B 0.56 0.0035 < 0.0001   2   2  23
    ##  ASBR08C ~/~ ASBR08C 0.90 0.0057 < 0.0001   3   3  24
    ##  ASBR08D ~/~ ASBR08D 0.71 0.0045 < 0.0001   4   4  25
    ##  ASBR08E ~/~ ASBR08E 0.71 0.0045 < 0.0001   5   5  26
    ##  ASBR08F ~/~ ASBR08F 0.70 0.0044 < 0.0001   6   6  27

``` r
# Look at the model fit
network3 %>% psychonetrics::fit()
```

    ##            Measure     Value
    ##               logl -84797.63
    ##  unrestricted.logl -84797.63
    ##      baseline.logl -97894.55
    ##               nvar         6
    ##               nobs        27
    ##               npar        27
    ##                 df       ~ 0
    ##          objective      2.59
    ##              chisq       ~ 0
    ##             pvalue         1
    ##     baseline.chisq  26193.85
    ##      baseline.npar        12
    ##        baseline.df        15
    ##    baseline.pvalue       ~ 0
    ##                nfi         1
    ##               pnfi       ~ 0
    ##                tli          
    ##               nnfi         1
    ##                rfi          
    ##                ifi         1
    ##                rni         1
    ##                cfi         1
    ##              rmsea          
    ##     rmsea.ci.lower       ~ 0
    ##     rmsea.ci.upper       ~ 0
    ##       rmsea.pvalue       ~ 0
    ##             aic.ll 169649.26
    ##            aic.ll2 169649.38
    ##              aic.x       ~ 0
    ##             aic.x2        54
    ##                bic 169849.87
    ##               bic2 169764.07
    ##            ebic.25 169898.25
    ##             ebic.5 169946.63
    ##            ebic.75 169985.33
    ##              ebic1 170043.38

``` r
# Compare the models
comparison <- psychonetrics::compare(
  `1. Original model`  = network2,
  `2. Sparse Model: Only Pruning` = network3)

print(comparison)
```

    ##                          model DF       AIC       BIC RMSEA Chisq Chisq_diff
    ##              1. Original model  0 169649.26 169849.87         ~ 0           
    ##  2. Sparse Model: Only Pruning  0 169649.26 169849.87         ~ 0        ~ 0
    ##  DF_diff p_value
    ##                 
    ##        0       1
    ## 
    ## Note: Chi-square difference test assumes models are nested.

# Exploratory Graph Analysis

``` r
# Dimension stability analysis via EGAnet
bootEGA1 <- EGAnet::bootEGA(
  data = scr, 
  cor = "cor_auto",
  uni.method = "louvain",
  iter = 500, # Number of replica samples to generate
  # resampling" for n random subsamples of the original data
  # parametric" for n synthetic samples from multivariate normal dist.
  type = "parametric", 
  # EGA Uses standard exploratory graph analysis
  # EGA.fit Uses total entropy fit index (tefi) to determine best fit of EGA
  # hierEGA Uses hierarchical exploratory graph analysis
  EGA.type = "EGA", 
  model = "glasso", 
  algorithm = "walktrap", # or "louvain" (better for unidimensional structures)
  # use "highest_modularity", "most_common", or "lowest_tefi"
  consensus.method = "highest_modularity", 
  typicalStructure = TRUE, # typical network of partial correlations
  plot.typicalStructure = TRUE, # returns a plot of the typical network
  ncores = 8, # Number of cores to use in computing results
  seed = 2024 # set the seed for replicability
)
```

![](network_files/figure-gfm/ch7-1.png)<!-- -->![](network_files/figure-gfm/ch7-2.png)<!-- -->

``` r
# View the number of communities
bootEGA1$EGA
```

    ## Model: GLASSO (EBIC with gamma = 0.5)
    ## Correlations: cor_auto
    ## Lambda: 0.0744050586187611 (n = 100, ratio = 0.1)
    ## 
    ## Number of nodes: 6
    ## Number of edges: 13
    ## Edge density: 0.867
    ## 
    ## Non-zero edge weights: 
    ##      M    SD    Min   Max
    ##  0.112 0.229 -0.180 0.540
    ## 
    ## ----
    ## 
    ## Algorithm:  Louvain
    ## 
    ## Number of communities:  1
    ## 
    ## ASBR08A ASBR08B ASBR08C ASBR08D ASBR08E ASBR08F 
    ##       1       1       1       1       1       1 
    ## 
    ## ----
    ## 
    ## Unidimensional Method: Louvain
    ## Unidimensional: Yes
    ## 
    ## ----
    ## 
    ## TEFI: 0

``` r
bootEGA1$typicalGraph$typical.dim.variables
```

    ##     items dimension
    ## 1 ASBR08A         1
    ## 2 ASBR08B         1
    ## 3 ASBR08C         1
    ## 4 ASBR08D         1
    ## 5 ASBR08E         1
    ## 6 ASBR08F         1

``` r
# Dimension (i.e., structural) stability results
dim_scr <- EGAnet::dimensionStability(bootEGA1)
```

![](network_files/figure-gfm/ch7-3.png)<!-- -->

``` r
dim_scr$dimension.stability
```

    ## $structural.consistency
    ## 1 
    ## 1 
    ## 
    ## $average.item.stability
    ## 1 
    ## 1

``` r
# Item stability results
dim_scr$item.stability
```

    ## EGA Type: EGA 
    ## Bootstrap Samples: 500 (Parametric)
    ## 
    ## Proportion Replicated in Dimensions:
    ## 
    ## ASBR08A ASBR08B ASBR08C ASBR08D ASBR08E ASBR08F 
    ##       1       1       1       1       1       1

``` r
dim_scr$item.stability$plot # to see only the plot
```

![](network_files/figure-gfm/ch7-4.png)<!-- -->

# Random-Intercept EGA

Random-Intercept EGA estimates the number of dimensions after
controlling for wording effects. EGA is applied to a residual
correlation matrix after subtracting a random intercept factor model
with equal unstandardized loadings from all the regular and unrecoded
reversed items in the data.

``` r
riEGA <- EGAnet::bootEGA(
  scr, 
  uni.method = "LE",
  iter = 500,
  type = "parametric",
  cor = "cor_auto",
  model = "glasso",
  EGA.type = "riEGA", # select random-intercept EGA here
  consensus.method  = "highest_modularity", 
  consensus.iter = 100,
  algorithm="walktrap",
  seed = 2024
)
```

![](network_files/figure-gfm/ch8-1.png)<!-- -->

``` r
print(riEGA)
```

    ## Model: GLASSO (EBIC)
    ## Correlations: cor_auto
    ## Algorithm:  Leading Eigenvector
    ## Unidimensional Method:  Leading Eigenvector
    ## 
    ## ----
    ## 
    ## EGA Type: riEGA 
    ## Bootstrap Samples: 500 (Parametric)
    ##              
    ##             1
    ## Frequency:  1
    ## 
    ## Median dimensions: 1 [1, 1] 95% CI

# Network Tree Analysis

In the final step, we will analyze whether the network structure of the
SCR scale differs by gender and reading achievement (i.e., above or
below the mean reading achievement score) using network tree analysis.

``` r
scr2 <- dplyr::select(scr_canada, starts_with("ASBR08"), IDCNTRY, ITSEX, ASRREA01) %>%
  dplyr::filter_all(all_vars(!is.na(.))) %>%
  dplyr::mutate_all(~na_if(., 9)) %>%
  dplyr::mutate(gender = as.factor(ifelse(ITSEX==1,"girl", 
                                          ifelse(scr_canada$ITSEX==2, "boy", NA))),
                reading = as.factor(ifelse(ASRREA01 < 500, "below average", "above average")),
                province = as.factor(ifelse(IDCNTRY == 9130, "Newfoundland and Labrador", ifelse(IDCNTRY == 9133, "Quebec", ifelse(IDCNTRY == 9134, "Alberta", "British Columbia"))))) %>%
  as.data.frame()


nmt_scr <- networktree::networktree(nodevars=scr, 
                                    splitvars=scr2[,c("gender", "reading", "province")],
                                    transform="pcor",
                                    na.action=na.omit)

plot_scr <- plot(nmt_scr, transform="pcor", maximum=0.2, edge.width=5, vsize=10, theme="colorblind", tnex = 3, partyargs=list(ep_args = list(justmin = 15), gp = grid::gpar(cex = .5)))
```

![](network_files/figure-gfm/ch9-1.png)<!-- -->
