# Dynamical Dark Energy and the Results of the DESI Collaboration

*Bachelor's Thesis in Physics · Universitat de València · 2025–2026*

Code and data analysis for my Bachelor's Thesis (*Trabajo de Fin de Grado*) in Physics at the Universitat de València.

* **Title:** *Dynamical Dark Energy and the Results of the DESI Collaboration* (*Energía Oscura Dinámica y los resultados de la Colaboración DESI*)
* **Author:** Pablo Ramón Muro
* **Supervisors:** Juan Herrero García and Maximilian Berbig (Departamento de Física Teórica, IFIC, CSIC–UV)
* **Degree:** Double Bachelor's Degree in Physics and Mathematics, Universitat de València
* **Academic year:** 2025–2026

> The written thesis is not included in this repository yet. It will be added once it has been formally submitted to the university.

## Overview

This repository contains the three Wolfram Language scripts used to produce every numerical result and every figure of the thesis:

1. the background evolution of a Friedmann–Lemaître–Robertson–Walker universe with different energy contents and a constant equation of state;
2. the $\chi^2$ fit of the $w_0w_a\mathrm{CDM}$ model (dynamical dark energy in the CPL parametrisation) to the DESI DR2 data, combined with BBN, a compressed CMB likelihood and compressed Pantheon+ supernovae;
3. the expansion history of the best-fit CPL cosmologies, compared with fiducial $\Lambda$CDM.

The scripts are self-contained: all the measurements that enter the fit are written explicitly inside the code, together with the reference of the paper they were taken from, so no external data files are needed.

## Requirements

* Wolfram Mathematica 12.0 or later. Only built-in functions are used (`NDSolveValue`, `NMinimize`, `FindMinimum`, `NIntegrate`, `ContourPlot`, `Export`), so no package has to be installed.
* No external datasets.

The scripts are stored as plain-text `.wl` files so that they can be read directly on GitHub. To run one of them, open it with Mathematica (`File > Open`) and evaluate the whole notebook (`Evaluation > Evaluate Notebook`). Each script calls `SetDirectory[NotebookDirectory[]]` before exporting, so the PDF figures are written to the same folder as the script.

## Repository structure

### `plots_DE_with_w_const.wl`

Background cosmology with a constant equation of state. It defines $E(z)=H(z)/H_0$, the deceleration parameter $q(z)$ and the cosmological distances for a fiducial $\Lambda$CDM model ($\Omega_M=0.31$, $\Omega_r=9\times10^{-5}$, $h=0.7$), and compares it with three limiting universes: Einstein–de Sitter, open and de Sitter.

It prints the characteristic numbers of the fiducial model (the matter–radiation and matter–$\Lambda$ equalities, the redshift at which the expansion starts to accelerate, $q_0$ and the age of the universe) and exports six figures:

| File | Content |
| --- | --- |
| `1_Historia_Expansion_Ez.pdf` | $E(z)=H/H_0$ for the four universes. |
| `2_Densidades_Eras.pdf` | Densities $\rho_i/\rho_{c,0}$ against the scale factor (log–log), with the radiation–matter and matter–$\Lambda$ equalities marked. |
| `3_Fracciones_Energia_Omega.pdf` | Energy fractions $\Omega_i(z)$ against $1+z$, with the same two equalities. |
| `4_Factor_Escala_at.pdf` | Scale factor $a(t)$ for the four universes, with $t-t_0=0$ at the present time. |
| `5_Desaceleracion_qz.pdf` | Deceleration parameter $q(z)$, with the onset of acceleration of $\Lambda$CDM marked. |
| `6_Distancias_Cosmologicas.pdf` | Comoving, Hubble, angular-diameter and luminosity distances in $\Lambda$CDM (logarithmic vertical axis), showing the maximum of $D_A$ and the duality $D_L=(1+z)^2D_A$. |

### `DESI_fit.wl`

Fit of the $w_0w_a\mathrm{CDM}$ model (CPL parametrisation $w(a)=w_0+w_a(1-a)$) to the DESI DR2 data in four combinations: BAO alone, BAO + BBN, BAO + CMB (compressed) and BAO + CMB + SNe (compressed Pantheon+). For each combination the script obtains:

* the best fit, that is, the minimum of $\chi^2$ and the parameters that reach it;
* the confidence region in the $(w_0,w_a)$ plane, built by profiling the remaining parameters;
* the position of $\Lambda$CDM inside that region and its exact confidence level (in sigmas and in % CL);
* the significance of $w_0w_a\mathrm{CDM}$ against $\Lambda$CDM, $\Delta\chi^2_{\mathrm{MAP}}\to\sigma$;
* the $1\sigma$ uncertainty of every parameter, obtained from the curvature of its one-dimensional profile.

It prints three tables (the parameters of the $\Lambda$CDM fits, the parameters of the $w_0w_a\mathrm{CDM}$ fits and the significances, all of them compared with the published DESI values) and exports seven figures:

| File | Content |
| --- | --- |
| `A_Contornos_BAO.pdf` | 68% and 95% confidence regions in the $(w_0,w_a)$ plane for BAO alone. |
| `B_Contornos_BAO_BBN.pdf` | The same regions for BAO + BBN. |
| `C_Contornos_BAO_CMB.pdf` | 68%, 95% and 99.7% regions for BAO + CMB. |
| `D_Contornos_BAO_CMB_SNe.pdf` | The same three regions for BAO + CMB + SNe. |
| `E_Tension_Hubble_H0.pdf` | $H_0$ inferred from the six fits (three datasets $\times$ two models), with SH0ES and Planck as external anchors. |
| `F_Restriccion_OmegaM_Verosimilitud.pdf` | One-dimensional likelihood of $\Omega_M$ in $\Lambda$CDM for the three datasets. |
| `G_Contornos_Solape_w0wa.pdf` | The BAO, BAO + CMB and BAO + CMB + SNe regions on the same axes. |

### `plots_CPL.wl`

Expansion history of the two best-fit CPL cosmologies obtained in `DESI_fit.wl`, compared with fiducial $\Lambda$CDM. The best-fit parameters are isolated in a single block at the top of the file, so they can be updated if the fit changes.

It prints the age of the universe, the redshift at which the expansion starts to accelerate and $q_0$ for each of the three models, and exports seven figures:

| File | Content |
| --- | --- |
| `A_Historia_Expansion_Ez.pdf` | $E(z)$ for the three cosmologies. |
| `B_Desaceleracion_qz.pdf` | $q(z)$ for the three cosmologies, with the $q=0$ line. |
| `C_Factor_Escala_at.pdf` | $a(t)$ for the three cosmologies, each one with its own $H_0$. |
| `D_Densidades_Eras_BAO_CMB.pdf`, `D_Densidades_Eras_BAO_CMB_SNe.pdf` | Densities $\rho_i/\rho_{c,0}$ of radiation, matter and dark energy for each best fit, with the two equalities marked. |
| `E_Fracciones_Energia_BAO_CMB.pdf`, `E_Fracciones_Energia_BAO_CMB_SNe.pdf` | Energy fractions $\Omega_i(z)$ for each best fit. |

## Data used in the fit

| Dataset | Source |
| --- | --- |
| BAO: 13 measurements of $D_M/r_d$, $D_H/r_d$ and $D_V/r_d$ with their covariance | DESI DR2, Table IV of [arXiv:2503.14738](https://arxiv.org/abs/2503.14738) |
| BBN: Gaussian prior on $\omega_b$ | Schöneberg et al. (2021), BBN + primordial deuterium |
| CMB: compressed distance priors $(R,\ell_A,\omega_b)$ with their full covariance | Chen, Huang & Wang, JCAP (2019), derived from Planck 2018 TT,TE,EE+lowE; $r_{\mathrm{drag}}$ and $r_*$ from Planck 2018 VI (Aghanim et al. 2020), Table 2 |
| SNe Ia: 7 compressed measurements of $1/E(z)$ with their $7\times7$ covariance | Table II of [arXiv:2408.17318](https://arxiv.org/abs/2408.17318), a compression of the full Pantheon+ catalogue |
| Reference values used for comparison | Eqs. (17), (19) and (24) and Tables 5 and VI of DESI DR2 |

The fitted values are compared throughout with the ones published by DESI. The CMB likelihood used here is the compressed one (three numbers), not the full Planck likelihood, so the results are close to the published ones but not identical.

## Notes

* The comments inside the code are written in English, but the axis labels, legends and messages printed by the scripts are in Spanish, because the figures are the ones included in the thesis, which is written in Spanish.
* The dark energy density in the CPL models is extrapolated outside the redshift range covered by the data ($z\gtrsim2$); the corresponding part of the curves should be read with that in mind.

## License

The code in this repository is released under the MIT License (see [`LICENSE`](LICENSE)).

## Contact

Pablo Ramón Muro — [pramonmuro@gmail.com](mailto:pramonmuro@gmail.com) · [LinkedIn](https://linkedin.com/in/pramonmuro)

---

## Resumen (castellano)

Este repositorio contiene los tres códigos de Wolfram Language utilizados en mi Trabajo de Fin de Grado en Física en la Universitat de València, *Energía Oscura Dinámica y los resultados de la Colaboración DESI*, dirigido por Juan Herrero García y Maximilian Berbig (IFIC, CSIC–UV). Con ellos se realiza el ajuste $\chi^2$ del modelo $w_0w_a\mathrm{CDM}$ a los datos públicos de DESI DR2 (BAO, BAO + BBN, BAO + CMB comprimido y BAO + CMB + SNe), se obtienen las regiones de confianza en el plano $(w_0,w_a)$ y la significancia frente a $\Lambda$CDM, y se generan todas las figuras de la memoria. Cada archivo comienza con un resumen de lo que hace y los comentarios explican el origen físico y numérico de cada paso.
