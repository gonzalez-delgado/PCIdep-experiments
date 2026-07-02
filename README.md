# PCIdep-experiments
[![SWH](https://archive.softwareheritage.org/badge/origin/https://github.com/gonzalez-delgado/PCIdep-experiments/)](https://archive.softwareheritage.org/browse/origin/?origin_url=https://github.com/gonzalez-delgado/PCIdep-experiments)

This code reproduces the numerical analyses and real data experiment in "Post-clustering Inference under Dependence", González-Delgado, Deronzier, Cortés and Neuvial 2023, [arxiv:2310.11822](https://arxiv.org/abs/2310.11822). Analyses make use of the $\texttt{R}$ package [PCIdep](https://github.com/gonzalez-delgado/PCIdep). This repository is organized as follows:

* The file [figure.1.R](https://github.com/gonzalez-delgado/PCIdep-experiments/blob/main/figure.1.R) reproduces the first numerical analysis described in the introduction of the paper, generating Figure 1.
* The file [figure.2.R](https://github.com/gonzalez-delgado/PCIdep-experiments/blob/main/figure.2.R) reproduces the second numerical analysis described in the introduction of the paper, generating Figure 2.
* The file [global.null.R](https://github.com/gonzalez-delgado/PCIdep-experiments/blob/main/global.null.R) reproduces the numerical analysis described in Section 4.1 and produces Figures 3 and G.2. This simulates the distribution of $p$-values under a global null hypothesis in various settings. 
* The file [sigma.overestimation.R](https://github.com/gonzalez-delgado/PCIdep-experiments/blob/main/sigma.overestimation.R) reproduces the numerical analyses described in Section 4.2, producing Figure 4 and Figure G.3. This illustrates that $p$-values are asymptotically super-uniform under the null hypothesis when one of the covariance matrices is asymptotically over-estimated.
* The file [conditional.power.R](https://github.com/gonzalez-delgado/PCIdep-experiments/blob/main/conditional.power.R) reproduces the conditional power analyses in the first part of Appendix C, producing Figure C.1.
* The file [conditional.power.spherical.covariance.R](https://github.com/gonzalez-delgado/PCIdep-experiments/blob/main/conditional.power.spherical.covariance.R) reproduces the conditional power analyses in the second part of Appendix C.1, producing Figure C.2.
* The file [non.CS.U.global.null.R](https://github.com/gonzalez-delgado/PCIdep-experiments/blob/main/non.CS.U.global.null.R) produces the first part of the analyses presented in Appendix D.1, producing Figures D.1 and G.4.
* The file [non.CS.U.overestimation.R](https://github.com/gonzalez-delgado/PCIdep-experiments/blob/main/non.CS.U.overestimation.R) reproduce the analyses in the second part of Appendix D.1, producing Figures D.2 and G.5.
* The file [non.admissible.U.overestimation.R](https://github.com/gonzalez-delgado/PCIdep-experiments/blob/main/non.admissible.U.overestimation.R) reproduces the analyses in Section D.2, producing Figures D.3 and G.6.
* The file [ignore.dependency.R](https://github.com/gonzalez-delgado/PCIdep-experiments/blob/main/ignore.dependency.R) reproduces the analyses in Section D.3, producing Figures D.4 and G.7.
* The file [matrix.t.perturbation.R](https://github.com/gonzalez-delgado/PCIdep-experiments/blob/main/matrix.t.perturbation.R) reproduces the analyses in Section D.4, producing Figures D.5 and G.8.
* The file [protein.clustering.R](https://github.com/gonzalez-delgado/PCIdep-experiments/blob/main/protein.clustering.R) reproduces experiment on real data of protein structures, presented in Section 5. The file
produces Table 1, Figure 5 and Figure G.8. The dataset can be downloaded [here](https://zenodo.org/doi/10.5281/zenodo.10021201). 
* The file [pgamma.sim.R](https://github.com/gonzalez-delgado/PCIdep-experiments/blob/main/pgamma.sim.R) reproduces the analysis of Section F.3, producing Figure F.1.
