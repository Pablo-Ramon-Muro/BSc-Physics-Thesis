# BSc-Physics-Thesis

Welcome! 👋 This repository contains all the code, notebooks, and datasets used for my Physics Bachelor's Thesis: *Dynamical Dark Energy and the results of the DESI Collaboration*. 

Here you will find the scripts used to perform the $\chi^2$ statistical analysis, the cosmological distance calculations, and the generation of the plots included in the thesis.

## Code structure / Estructura del código

* `plots_DE_with_w_const.wl`: Mathematica notebook used for the calculations and plots of the LCDM model with different energy content.
* `DESI_fit.wl`: Mathematica notebook where we perform a fit of the $w_0 w_a \text{CDM}$ model (dynamical dark energy, CPL) to the DESI DR2 data in four combinations: BAO alone, BAO + BBN, BAO + CMB (compressed), and BAO + CMB + SNe (compressed Pantheon+). For each one we obtain:
    - The best fit (minimum of $\chi^2$).
    - The confidence region in the $(w_0, w_a)$ plane.
    - The position of LCDM and its exact confidence value (sigmas, % CL).
    - The significance of $w_0 w_a \text{CDM}$ vs LCDM ($\Delta\chi^2_{\text{MAP}} \to$ sigmas).
    - Comparison tables with DESI, all framed and in the usual format.
* `plots_CPL.wl`: Mathematica notebook where we obtain plots of the expansion history comparing three cosmologies, using the CPL parametrisation $w(a) = w_0 + w_a(1-a)$:
    - Fiducial LCDM ($\Omega_M=0.31$, $w=-1$).
    - $w_0 w_a \text{CDM}$ BAO+CMB (best fit).
    - $w_0 w_a \text{CDM}$ BAO+CMB+SNe (best fit).
