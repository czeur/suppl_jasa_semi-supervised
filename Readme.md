# Guidance for the codes

This file provides a guidance for the codes used to produce the empirical results in the paper **Extreme value statistics in semi-supervised models** by *Ahmed, Einmahl and Zhou*. The codes contain two parts: simulation and application. For the application part, additional data files are needed, which are not provided together with the codes. The data are available upon request, with the consent of Meteo France.

All codes are written in R. They can be run under R 4.3.0 (2023-04-21).

## Simulation

### Figure 1 with varying k

The codes to create plots in Figure 1 are given in the file `Figure 1 code.R`. Running this codes leads to the desired Figure.

### Variance reduction with pre-selected k
The codes to produce simulation results for the variance reduction are in the files `Paper_code_JASA.R` and `Paper_code_g_t=0_JASA.R`. The difference is that the first file handles non-zero extreme value index, while the second file handles gamma=0 specifically. For the rest, the two files have the same structure.

In the input part, one can specify the following inputs 

```
repl=500
results_tab=matrix(NA,repl,9)
n=500
m=1000
k=125
p=(1/n) 
g_t=-0.3 #true value of gamma
```
In addition, in the part "Generating simulation data", the remark guides how to specify different data generating process.

The results with different model specifications and parameter choices are saved in `result.xlsx`. The results include both the 2d and 3d results for the variance reduction. These results are required for generating Table 1, Figure 2, Figure 3, Table 3 and Table 4.

### Additional simulation results in supplementary material

The codes to produce simulation results from the biavariate Gaussian copula in the supplementary material are in the file `Paper_code_g_t=0_JASA_normal.R`. 

Running this file will provide simulation results for both 2d and 3d normal distributions. The results are saved in `result.xlsx`. Only the 2d results are reported in the supplementary material.

## Application

The application codes generate Table 5 and Figure 4 in the paper. The file `semisupervise.R` contains relevant functions. Running the file `application.R` leads to the desired Table and Figure. 

Running this file requires additional data files that are not included in the codes. They are available upon request with the consent of Meteo France.

Please check the remarks in those two files for further information. Part of the remarked codes can help to produce a Hill plot with varying k. However, to produce the final Table and Figure, we need to fix k as its current level.

## Simulations