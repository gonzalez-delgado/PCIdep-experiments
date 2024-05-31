# PCIdep-experiments
This code reproduces the numerical analyses and real data experiment in "Post-clustering Inference under Dependency", González-Delgado, Cortés and Neuvial 2023, [arxiv:2310.11822](https://arxiv.org/abs/2310.11822). Analyses make use of the $\texttt{R}$ package [PCIdep](https://github.com/gonzalez-delgado/PCIdep). This repository is organized as follows:

* The file [figure.1.R](https://github.com/gonzalez-delgado/PCIdep-experiments/blob/main/figure.1.R) reproduces the first numerical analysis described in the introduction of the paper, generating Figure 1.
* The file [figure.2.R](https://github.com/gonzalez-delgado/PCIdep-experiments/blob/main/figure.2.R) reproduces the second numerical analysis described in the introduction of the paper, generating Figure 2.
* The file [global.null.R](https://github.com/gonzalez-delgado/PCIdep-experiments/blob/main/global.null.R) reproduces the numerical analysis described in Section 5.1. This simulates the distribution of $p$-values under a global null hypothesis in various settings. The script generates Figures 3 and 10 in the paper.
* The file [sigma.overestimation.R](https://github.com/gonzalez-delgado/PCIdep-experiments/blob/main/sigma.overestimation.R) reproduces the numerical analyses described in Section 5.2. This illustrates that $p$-values are asymptotically super-uniform under the null hypothesis when one of the covariance matrices is asymptotically over-estimated. The script produces Figures 4 and 13 in the paper.
* The file [conditional.power.R](https://github.com/gonzalez-delgado/PCIdep-experiments/blob/main/conditional.power.R) reproduces the conditional power analyses in Section 5.3, when both covariance matrices are known and when one of them is asymptotically over-estimated. This generates Figure 5.
* The file [non.admissible.U.R](https://github.com/gonzalez-delgado/PCIdep-experiments/blob/main/non.admissible.U.R) reproduces the analyses in Section 5.4.1, generating Figures 6 and 14.
* The file [ignore.dependency.R](https://github.com/gonzalez-delgado/PCIdep-experiments/blob/main/ignore.dependency.R) reproduces the analyses in Section 5.4.2, generating Figures 7 and 15.
* The file [protein.clustering.R](https://github.com/gonzalez-delgado/PCIdep-experiments/blob/main/protein.clustering.R) reproduces experiment on real data of protein structures, presented in Section 6. The dataset can be downloaded [here](https://zenodo.org/doi/10.5281/zenodo.10021201). This script produces Figures 8 and 9.
* The file file [section.D2.R](https://github.com/gonzalez-delgado/PCIdep-experiments/blob/main/section.D2.R) reproduces the analysis described in Section D.3, generating Figures 11 and 12. 
