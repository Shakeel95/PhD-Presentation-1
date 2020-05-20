# PhD-Presentation-1

This repository contains code for building my my first year PhD presentation.

## Abstract 

In many settings prominent features of univariate time series can be described  well  by a stylised model in which weakly dependent noise fluctuates around a piecewise linear trend. When experimental setups lead to many such time series being collected together as a panel, subgroups will often exhibit visual similarity. To motivate high dimensional trend segmentation, the first part of the talk will introduce a data example in which classical time series clustering methods fail the "visual similarity test" while a simple procedure based on changepoints performs surprisingly well.

The second part of the talk will review existing methods for univariate trend segmentation, and discuss strategies for extending to the multivariate setting. While the problem of high dimensional trend segmentation is relatively new, the problem of identifying prominent features or *landmarks* in a sample of curves dates back to o Kneip & Gasser (1992) . In functional data analysis landmarks are used to perform curve registration. In the context of high dimensional time series, feature misalignment may in fact be informative as it can be interpreted as reflecting  components of a connected system responding to impulses at different speeds; analogous to dynamic factor models. 

The final part of the talk will propose a factor model framework for high dimensional trend segmentation problems. We expect such a model to exhibit a *blessing of dimensionality*, in the sense that estimates for the locations of aligned changepoints improve with T and n. Estimating the model amounts to recovering a common basis for piecewise linear functions on [1,T]. This will be compared to the associated eigen-basis problem which arises when estimating parameters in a standard factor model, as well as to the problem of principal component analysis for functional data. 


## Usage

Depends on R (>= 3.6). To make sure all packages needed to execute scrips in the repository are installed do: 

``` R 
str(ip <- installed.packages()) 
rp <- scan("required.txt", what = "", sep = "\n")

sapply(rp, function(i) if (!(i %in% rownames(ip))) install.packages(i))
```
