(* ::Package:: *)

(* =====================================================================
   DESI_fit.wl   (Wolfram Language / Mathematica)
   ---------------------------------------------------------------------
   Part of the code of the Bachelor's Thesis in Physics "Dynamical Dark
   Energy and the results of the DESI Collaboration" (Pablo Ram\[OAcute]n Muro,
   Universitat de Valencia).

   What this script does
   ---------------------
   It fits the w0waCDM model (dynamical dark energy with the CPL
   parametrisation w(a) = w0 + wa (1-a)) to the DESI DR2 data in four
   combinations:
        (1) BAO alone,      (2) BAO + BBN,
        (3) BAO + CMB (compressed),
        (4) BAO + CMB + SNe (compressed Pantheon+).
   For each combination it obtains:
        - the best fit, that is, the minimum of chi^2 and the parameters
          that reach it,
        - the confidence region in the (w0,wa) plane,
        - the position of LCDM inside that region and its exact confidence
          level (in sigmas and in % CL),
        - the significance of w0waCDM against LCDM (Dchi2_MAP -> sigmas),
        - the 1 sigma uncertainty of every parameter,
   and it builds the comparison tables with the values published by DESI.

   Data
   ----
   Everything that enters the chi^2 is written explicitly in Block 0, with
   the reference of the paper it was taken from, so the script needs no
   external file: the 13 BAO measurements of DESI DR2 with their covariance,
   the BBN prior on omega_b, the compressed CMB distance priors of Planck
   2018 and the compressed Pantheon+ supernovae.

   Structure of the file
   ---------------------
   Block 0  : constants and data.
   Block 1  : the model E(z) and the comoving integral.
   Block 2  : the chi^2 of each combination.
   Block 3  : best fit and significance of each case.
   Block 3b : common format of the figures.
   Block 4  : confidence regions in the (w0,wa) plane.
   Block 4b : H0 and the Hubble tension.
   Block 4c : constraint on Omega_M from the one-dimensional likelihood.
   Block 4d : the three confidence regions on the same axes.
   Block 5  : comparison tables with DESI.
   Block 6  : summary of the comparison with the published results.
   Block 7  : how the supernovae are included.
   Block 8  : export of the figures to PDF.

   Output
   ------
   Three tables printed in the notebook and seven figures exported to
   vector PDF in the folder of the notebook.

   Mathematica functions that do the work
   --------------------------------------
   NMinimize    finds the global minimum inside a box of parameters.
   FindMinimum  finds the local minimum starting from a given point.
   NDSolveValue solves the comoving integral as an ODE (see Block 1).
   NIntegrate   integrates numerically.
   ContourPlot  draws the contour lines of a grid of numbers.
   Grid[..., Frame -> All]  builds the framed tables.
   ===================================================================== *)


(* =====================================================================
   Block 0 : constants and data from Table IV (DESI DR2, 2503.14738)
   ===================================================================== *)
cc  = 299792.458;         (* speed of light [km/s]. *)
wg  = 2.47*^-5;           (* omega_gamma = Omega_gamma h^2  (photons).      *)
wr  = 4.178*^-5;          (* omega_r     = Omega_r     h^2  (photons+nu).   *)
wnu = 0.0006;             (* omega_nu    (Sum m_nu = 0.06 eV, NH minimum).  *)

(* --- Data vector d ---
   Each DESI bin contributes either one isotropic measurement DV/rd (low signal),
   or two anisotropic measurements (DM/rd , DH/rd). We stack them all in a
   single vector d, and store in parallel the effective redshift (zvec) and
   the type of each component (tvec: "M"=DM/rd, "H"=DH/rd, "V"=DV/rd).    *)
dvec = {7.942,                    (* BGS        DV/rd  (isotropic).         *)
        13.588, 21.863,           (* LRG1       DM/rd , DH/rd .             *)
        17.351, 19.455,           (* LRG2       DM/rd , DH/rd .                        *)
        21.576, 17.641,           (* LRG3+ELG1  DM/rd , DH/rd .                           *)
        27.601, 14.176,           (* ELG2       DM/rd , DH/rd .                            *)
        30.512, 12.817,           (* QSO        DM/rd , DH/rd .                           *)
        38.988,  8.632};          (* Lya        DM/rd , DH/rd .                           *)
zvec = {0.295, 0.510,0.510, 0.706,0.706, 0.934,0.934,
        1.321,1.321, 1.484,1.484, 2.330,2.330};
tvec = {"V", "M","H", "M","H", "M","H", "M","H", "M","H", "M","H"};

(* --- Covariance matrix C ---
   DESI gives the data with their covariance. Within the same bin, DM/rd and
   DH/rd are correlated (coefficient r, column r_MH of the table):
   their 2x2 block is {{sM^2, r sM sH},{r sM sH, sH^2}}. Between different bins
   the correlation is negligible (shot-noise dilution, see the guide),
   so C is block-diagonal: a 1x1 for BGS and a 2x2 for each anisotropic
   bin. The chi^2 will use C^{-1}. We take the data from Table IV of DESI DR2. *)
cov2[sM_, sH_, r_] := {{sM^2, r sM sH}, {r sM sH, sH^2}};(* 2x2 covariance matrices of each bin, except BGS. *)
bloques = {{{0.075^2}},                          (* BGS (1x1 block).        *)
   cov2[0.167, 0.425, -0.459],   (* LRG1 *)      cov2[0.177, 0.330, -0.404], (* LRG2 *)
   cov2[0.152, 0.193, -0.416],   (* LRG3+ELG1 *) cov2[0.318, 0.221, -0.434], (* ELG2 *)
   cov2[0.760, 0.516, -0.500],   (* QSO  *)      cov2[0.531, 0.101, -0.431]};(* Lya  *) (* Blocks that will generate the covariance
    matrix of all bins, diagonal with the blocks of each bin. *)

(* Function that assembles the block-diagonal covariance (13x13) from the list. *)
blockDiag[bl_] := Module[{n = Total[Length /@ bl], m, i = 1}, (*Dimension of the nxn matrix, with n the sum of those of all the blocks. *)
   m = ConstantArray[0., {n, n}]; (* Generates an nxn matrix of zeros. *)
   Do[m[[i ;; i + Length[b] - 1, i ;; i + Length[b] - 1]] = b; i += Length[b], {b, bl}];
   m]; (* This is a loop, where b is the element of the loop we are at. i += Length[b] makes the new i be the previous one + Length[b].
   The loop places the blocks defined before where they belong, to build the covariance matrix of all the bins. *)
Cmat = blockDiag[bloques];      (* We apply the function defined above to build the covariance matrix of all the bins. *)
Cinv = Inverse[Cmat];           (* Its inverse, C^{-1}, for the chi^2.     *)

(* --- Compressed CMB data: 3x3 "distance priors" (Planck 2018) ---
   Vector v=(R, l_A, omega_b). Means, sigmas and correlation matrix from
   Chen, Huang & Wang 2019 (JCAP), derived from Planck 2018 TT,TE,EE+lowE.
   I use the wCDM row of the table, not the LCDM one: since here I fit
   dynamical dark energy (w0waCDM), it is more consistent to take the priors
   derived allowing w != -1 (wCDM, the most general row of the table) than
   those of LCDM (fixed w=-1). The effect is small, because the distance
   priors are almost model-independent (Chen-Huang-Wang): relative to LCDM
   only R, l_A, wb and rho(R,l_A) change (0.46 -> 0.47); the cross terms
   rho(R,wb)=-0.66 and rho(l_A,wb)=-0.34 and the sigmas stay the same.
   The covariance is C = S.rho.S (S=diag(sigmas)); its inverse enters chi2_CMB.
   Note the cross terms rho(R,wb)=-0.66 and rho(l_A,wb)=-0.34, which are not
   diagonal. *)
vObsCMB = {1.7493, 301.462, 0.02239};             (* Observed means (wCDM row): v=(R, l_A, omega_b). *)
sigCMB  = {0.0046, 0.090, 0.00015};               (* Sigmas of the observed measurements. *)
rhoCMB  = {{1., 0.47, -0.66},                     (* Correlation matrix rho_ij (wCDM). *)
           {0.47, 1., -0.34},
           {-0.66, -0.34, 1.}};
CinvCMB = Inverse[DiagonalMatrix[sigCMB] . rhoCMB . DiagonalMatrix[sigCMB]]; (*Inverse of the CMB covariance matrix.*)

(* --- Compressed SNe Ia data: Pantheon+ (Table II of 2408.17318v3) --- 
   The Pantheon+ catalogue (1701 SNe) is compressed to 7 measurements of 1/E(z)
   at 7 redshift nodes: 1/E(z) is modelled as a cubic spline fitted to the full
   Pantheon+ and the script-M nuisance (mix of M_B and H0) is marginalised. The
   result is a Gaussian in {1/E(z_i)}: means eSN, sigmas sSN and correlation
   matrix rhoSN (7x7). It gives a chi^2_SNe that only depends on (Om,w0,wa), not
   on H0, rd, wb or M. Nodes and values exactly as in Table II of the paper.  *)
zSN = {0.07, 0.20, 0.35, 0.55, 0.90, 1.50, 2.00};        (* Redshift of the 7 nodes.   *)
eSN = {0.947, 0.8982, 0.815, 0.664, 0.790, 0.18, 0.46};  (* Means of 1/E(z_i).         *)
sSN = {0.014, 0.0095, 0.020, 0.024, 0.066, 0.11, 0.22};  (* Sigmas of 1/E(z_i).        *)
rhoSN = {{ 1.00, 0.15, 0.59, 0.10, 0.15, -0.05, 0.07},   (* 7x7 correlation matrix.    *)
         { 0.15, 1.00, -0.16, 0.20, 0.04, 0.04, 0.03},
         { 0.59, -0.16, 1.00, -0.24, 0.21, -0.14, 0.12},
         { 0.10, 0.20, -0.24, 1.00, -0.50, 0.32, -0.11},
         { 0.15, 0.04, 0.21, -0.50, 1.00, -0.66, 0.25},
         { -0.05, 0.04, -0.14, 0.32, -0.66, 1.00, -0.34},
         { 0.07, 0.03, 0.12, -0.11, 0.25, -0.34, 1.00}};
CinvSN = Inverse[DiagonalMatrix[sSN] . rhoSN . DiagonalMatrix[sSN]]; (* Inverse of the 7x7 SNe covariance. *)


(* =====================================================================
   Block 1 : the model  E(z)=H(z)/H0  and the comoving integral
   ---------------------------------------------------------------------
   w0waCDM (flat) with CPL  w(a)=w0+wa(1-a):
     E(z)^2 = orr(1+z)^4 + Om(1+z)^3 + Ode (1+z)^{3(1+w0+wa)} e^{-3 wa z/(1+z)}
   with Ode = 1-Om-orr, and 'orr' = Omega_r = omega_r/h^2 (radiation fraction
   today).
   ===================================================================== *)
Ez[z_, Om_, w0_, wa_, orr_] :=
   Sqrt[orr (1+z)^4 + Om (1+z)^3
        + (1 - Om - orr) (1+z)^(3 (1 + w0 + wa)) Exp[-3 wa z/(1+z)]];

(* --- Comoving integral  I(z) = Integral_0^z dz'/E(z')  ---
   This is the key to the speed of the script. Instead of calling NIntegrate
   over and over (slow: for the CMB one would have to integrate up to z~1090
   tens of thousands of times, since one keeps trying different values of the
   free parameters in order to minimise), I solve the integral only once as
   an ODE with NDSolve:
        I'(z) = 1/E(z),  I(0)=0,
   and I get an interpolating function I(z) valid for all z at once
   (including z_star). Reusing it for the 13 BAO observables and for the
   D_M at z_star turns about 7 integrals per evaluation into a single ODE: it
   is between 20 and 50 times faster and gives the same number.
   'zmax' = how far I need I(z): 2.5 for BAO, about 1100 if there is CMB.   *)
Ifun[Om_?NumericQ, w0_?NumericQ, wa_?NumericQ, orr_?NumericQ, zmax_?NumericQ] :=
   NDSolveValue[{u'[z] == 1/Ez[z, Om, w0, wa, orr], u[0] == 0}, u, {z, 0, zmax}]; (*u(0) is 0 because it would be integrating from 0 to 0. *)

(* (Iof is kept in case a single one-off integral is needed; it is not used in
    the fast loops.) *)
Iof[z_?NumericQ, Om_?NumericQ, w0_?NumericQ, wa_?NumericQ, orr_?NumericQ] :=
   NIntegrate[1/Ez[zp, Om, w0, wa, orr], {zp, 0, z}];

(* --- "Scale-free model" vector g, from an already computed I(z) ---
   The BAO observables are  DM/rd = A*gM , DH/rd = A*gH , DV/rd = A*gV ,
   with the same scale constant  A = c/(100 h rd) (see Block 2), and
       gM(z) = I(z),   gH(z) = 1/E(z),   gV(z) = ( z I(z)^2 / E(z) )^{1/3}.
   'Iz' is the interpolating function returned by Ifun.                     *)
gvecOf[Iz_, Om_?NumericQ, w0_?NumericQ, wa_?NumericQ, orr_?NumericQ] :=
   MapThread[
     Which[
       #2 == "M", Iz[#1],
       #2 == "H", 1/Ez[#1, Om, w0, wa, orr],
       True,      (#1 (Iz[#1])^2 / Ez[#1, Om, w0, wa, orr])^(1/3)
     ] &, {zvec, tvec}]; (* Depending on the tvec label, for each component it computes DV/rd, DM/rd or DH/rd, using the
     corresponding zeff, which is given by the zvec component. The result is a vector: [DV/rd, (DM/rd, DH/rd)_i]/A, with i
     running over all the bins.*)


(* =====================================================================
   Block 2 : the chi^2 of the four cases
   ---------------------------------------------------------------------
   The chi^2 measures the data-model distance weighted by the covariance:
        chi^2(param) = (d - t(param))^T . C^{-1} . (d - t(param)).
   The best fit is the 'param' that minimises this number.
   ===================================================================== *)

(* --- (1) BAO alone ----------------------------------------------------
   Here rd and h are not known separately: BAO only measures the combination
   A = c/(100 h rd). Since the model t = A*g is linear in A (and g does not
   depend on A), the minimum in A is obtained in closed form, analytically
   (differentiate chi^2 with respect to A and set it to 0), without affecting
   the value of the other parameters:
        A_opt = (g^T C^{-1} d)/(g^T C^{-1} g),
        chi^2_min(A) = d^T C^{-1} d - (g^T C^{-1} d)^2/(g^T C^{-1} g).
   It is exactly the same idea as marginalising the M of the supernovae.
   So NMinimize only has to search in (Om,w0,wa): A comes for free.        *)
chi2BAO[Om_?NumericQ, w0_?NumericQ, wa_?NumericQ] :=
   Module[{orr = wr/0.49, Iz, g, gCg, gCd, dCd}, (* \[CapitalOmega]_r = \[Omega]_r/h\.b2 and although we do not know h a priori, since wrad has so little influence
   at the z we work with, we can give it a very wide range of values that does not change the result. We give it the mean value between two extreme
   values of the Hubble tension.*)
     Iz  = Ifun[Om, w0, wa, orr, 2.5];       (* One ODE up to z=2.5.         *)
     g   = gvecOf[Iz, Om, w0, wa, orr];
     gCg = g . Cinv . g;  gCd = g . Cinv . dvec;  dCd = dvec . Cinv . dvec;
     dCd - gCd^2/gCg];                       (*dvec was the vector of observed BAO data. $$\chi^2(A) =
     (\mathbf{d}^T \mathbf{C}^{-1} \mathbf{d}) - 2A(\mathbf{g}^T \mathbf{C}^{-1} \mathbf{d}) + A^2(\mathbf{g}^T \mathbf{C}^{-1} \mathbf{g})$$.
     I define: dCd $= \mathbf{d}^T \mathbf{C}^{-1} \mathbf{d}
     $gCd $= \mathbf{g}^T \mathbf{C}^{-1} \mathbf{d}
     $gCg $= \mathbf{g}^T \mathbf{C}^{-1} \mathbf{g}$.
      A_{opt} = \frac{\text{gCd}}{\text{gCg}}.
      chi^2 already minimised in A: $$\chi_{min}^2BAO = \text{dCd} - \frac{(\text{gCd})^2}{\text{gCg}}$$. *)
Aopt[Om_?NumericQ, w0_?NumericQ, wa_?NumericQ] :=
   Module[{orr = wr/0.49, Iz, g},
     Iz = Ifun[Om, w0, wa, orr, 2.5]; g = gvecOf[Iz, Om, w0, wa, orr];
     (g . Cinv . dvec)/(g . Cinv . g)];                                 (*although we had already put them into the chi^2 formula, to store 
     the value of the
      optimal A in case we want it later, to get H_0*r_d, we compute it with this function.*)

(* --- rd from the DESI fitting formula (their eq. 2), valid only with BBN or
   CMB, where wb is a free parameter ------------------------------------
   Instead of integrating the sound horizon, DESI gives a calibrated fit:
        rd ~ 147.05 (wb/0.02236)^-0.13 (wbc/0.1432)^-0.23 (Neff/3.04)^-0.1 [Mpc]
   with wb=Omega_b h^2, wbc=Omega_cb h^2 (baryons+CDM), Neff=3.044.        *)
rdForm[wb_, wbc_, Neff_:3.044] :=
   147.05 (wb/0.02236)^-0.13 (wbc/0.1432)^-0.23 (Neff/3.04)^-0.1;

(* --- (2) BAO+BBN  and  (3) BAO+CMB ------------------------------------
   When adding BBN or CMB, rd becomes calibrated through (wb,wbc) and then A is
   no longer free: A = c/(100 h rd), with h a real parameter. We cannot minimise
   A without taking the rest into account, because A depends on rd, which in turn
   depends on wb, which is present there but also in the chi^2 of BBN or CMB. h
   cannot be minimised either, because it is present in the calculation
   wcb = Omega_M*h^2 - wnu, from which r_d is obtained. The base chi^2 is the BAO
   one with that A, plus the information from the extra data:
     - BBN: wb is a free parameter (baryon density) that calibrates rd. There is
            no extra term: the BBN information enters entirely through rd(wb).
     - Compressed CMB: I use the 3x3 distance priors of Planck 2018 with their
            full covariance (including the cross terms), the rigorous version of
            Chen, Huang & Wang 2019 (JCAP). The vector is v = (R, l_A, omega_b):
              R   = sqrt(Om) * I(z_star)          (shift parameter)
              l_A = pi * D_M(z_star) / r_star     (acoustic scale)
              omega_b (free parameter)
            with I(z_star)=Integral_0^z* dz/E, D_M(z_star)=(c/100h) I(z_star),
            r_star = 0.9819 rd (calibrated to centre l_A) and z_star~1090.
            Observed values and covariance (Planck 2018, Chen-Huang-Wang 2019;
            wCDM row, see Block 0):
              v_obs = (1.7493, 301.462, 0.02239)
              sigmas= (0.0046, 0.090,  0.00015)
              correlations: rho(R,l_A)=0.47, rho(R,wb)=-0.66, rho(l_A,wb)=-0.34
            chi2_CMB = (v_th - v_obs)^T . CinvCMB . (v_th - v_obs), see below.
            (Before I used (theta_*, wbc) diagonal; this replaces it with the
             full version. It gives the same, about 2.7 sigma, but it is the
             rigorous one.)

   How wb is treated: wb is now a free parameter (nuisance) with its prior:
     - BBN: Gaussian prior  (wb-0.02218)^2/0.00055^2  added to the chi^2
            (BBN + primordial deuterium D/H; Schoneberg et al. 2021).
     - CMB: no separate prior, because wb is the third component of the distance
            vector (mean 0.02239, sigma 0.00015; wCDM row) and its prior plus the
            cross terms with R and l_A are already inside the 3x3 C (Planck 2018).
   By leaving wb free, the uncertainty of rd(wb) propagates by itself to the error
   of H0, which is why H0 comes out with a more realistic error (for example with
   BBN about +-0.46 instead of +-0.23). The (w0,wa) contours barely change; the
   CMB significance drops a little (from about 2.7 to about 2.3 sigma) because
   marginalising wb relaxes the LCDM fit. The bound on wb in each NMinimize is
   only the search box.                                                     *)
chi2Full[Om_?NumericQ, h_?NumericQ, w0_?NumericQ, wa_?NumericQ, wb_?NumericQ, caso_] :=
   Module[{orr = wr/h^2, zmax, Iz, wbc, rd, A, g, res, x2, Istar, DMstar, Rth, lAth, dv},
     zmax = If[caso == "cmb", 1100., 2.5];  (* Only the CMB needs z~1090; otherwise we just go up to 2.5. *)
     Iz  = Ifun[Om, w0, wa, orr, zmax];     (* One ODE, reused for everything. *)
     wbc = Om h^2 - wnu;               (* Matter (b+CDM) minus neutrinos.      *)
     rd  = rdForm[wb, wbc];            (* Standard ruler calibrated with wb.   *)
     A   = cc/(100 h rd);              (* The scale is no longer free.         *)
     g   = gvecOf[Iz, Om, w0, wa, orr];
     res = dvec - A g; (* Observations minus model. *)
     x2  = res . Cinv . res;               (* BAO chi^2.                       *)
     If[caso == "bbn",                 (* BBN: Gaussian prior on wb.           *)
        x2 += (wb - 0.02218)^2/0.00055^2]; (*you add the bbn chi^2 to x2, which was the bao chi^2.*)
     If[caso == "cmb",                 (* CMB = 3x3 distance priors, see above. *)
        Istar  = Iz[1089.8];                             (* I(z_star)=Integral dz/E. *)
        DMstar = cc/(100 h) Istar;                       (* D_M(z_star).           *)
        Rth    = Sqrt[Om] Istar;                         (* R  = sqrt(Om)H_0*c I(z_star)/(H_0*c)= Sqrt[Om] Istar.*)
        (* r_star = (r*/rd) rd, with r*/rd = 144.43/147.09 = 0.9819 from the
           derived values of Planck 2018 VI (Aghanim et al. 2020), Table 2:
           r_drag=147.09 Mpc, r_*=144.43 Mpc.
           Justification: r_d and r_* are integrals from z_drag/z_star to infinity,
           that is, quantities of the early universe (z>1000), where dark energy
           is negligible. That is why the ratio does not depend on w0,wa (it is the
           same in LCDM as in w0waCDM). It only depends on wb,wcb, and I checked
           that it varies less than 0.02% over the whole fit range, so fixing it
           is almost exact. *)
        lAth   = Pi DMstar/(0.9819 rd);                  (* l_A = pi D_M/r_star.    *)
        dv     = {Rth, lAth, wb} - vObsCMB;              (* Residual (R,l_A,wb).    *)
        x2 += dv . CinvCMB . dv];                        (* chi^2 bao + chi^2 cmb, with the 3x3 quadratic form. *)
     x2]; (*that is the value that chi^2 returns to me.*)

(* --- (4) BAO+CMB+SNe: compressed SNe added to the CMB case --- 
   chi2SNe only depends on (Om,w0,wa): the model at each node is g_i=1/E(z_i) and
   it is compared with eSN through CinvSN. I use orr=wr/0.49 (the same as chi2BAO)
   because at z<=2 radiation is negligible, so chi2SNe does not need h.
   The combined chi^2, chi2CS, is that of BAO+CMB (chi2Full with "cmb", which
   already includes BAO+CMB with its wb) plus the SNe term; it shares the same
   5 parameters (Om,h,w0,wa,wb), the SNe do not add any new one.                *)
chi2SNe[Om_?NumericQ, w0_?NumericQ, wa_?NumericQ] := Module[{orr = wr/0.49, g},
   g = (1/Ez[#, Om, w0, wa, orr]) & /@ zSN;       (* Model 1/E(z_i) at the 7 nodes. *)
   (g - eSN) . CinvSN . (g - eSN)];               (* 7x7 SNe quadratic form.        *)
chi2CS[Om_?NumericQ, h_?NumericQ, w0_?NumericQ, wa_?NumericQ, wb_?NumericQ] :=
   chi2Full[Om, h, w0, wa, wb, "cmb"] + chi2SNe[Om, w0, wa];  (* BAO+CMB  +  SNe. *)

(* ---------------------------------------------------------------------
   From Dchi^2 to sigmas (Wilks' rule, 2 degrees of freedom = w0 and wa).
   A contour at Dchi^2 encloses a probability fraction
        CL(Dchi^2) = 1 - e^{-Dchi^2/2}          (chi^2 with 2 d.o.f.)
   and the equivalent Gaussian number of sigmas (two-tailed) is
        N(Dchi^2) = sqrt(2) * InverseErfc( e^{-Dchi^2/2} ).
   The argument 'd' must be positive (the Dchi^2 of a point above the
   minimum). I use Max[0,d] as a protection: if a fit does not converge and
   d<0 comes out (physically impossible), it returns 0 instead of nonsense,
   since InverseErfc only accepts [0,2]. A suspicious 0 warns that the fit
   failed. *)
clDe [d_] := 1 - Exp[-Max[0, d]/2];               (* d is the input value of the function, which will be Dchi_c^2. So,
we get prob(Dchi<= Dchi_c), which tells us in which confidence region a point lies       *)
sigDe[d_] := Sqrt[2] InverseErfc[Exp[-Max[0, d]/2]]; (* Equivalent sigmas. *)


(* =====================================================================
   Block 3 : best fit and significance of each case
   =====================================================================
   How I move through the parameters and get the minimum.
   NMinimize does not try random points without criterion: for a problem with
   bounds it uses a global optimiser (by default, differential evolution /
   Nelder-Mead). It keeps a population of candidates (Om,w0,wa,...),
   evaluates chi^2 at each one, combines and moves the best ones towards regions
   of lower chi^2, and repeats until convergence. I give it:
     - the function to minimise (the chi^2 of the case),
     - a search box (the flat priors: the range of each parameter),
     - a physical constraint: w0+wa<0, so that dark energy does not dominate in
       the early universe, which is the condition DESI imposes.
   It returns  {chi^2_minimum , {Om->..., w0->..., ...}}.

   Significance of w0waCDM against LCDM.
   LCDM is w0waCDM with (w0,wa)=(-1,0) fixed. I fit the two models
   separately and compare their minima:
        Dchi^2_MAP = chi^2_min(LCDM) - chi^2_min(w0waCDM)  >= 0.
   The larger it is, the more w0waCDM improves the fit. That Dchi^2 (2 d.o.f.)
   translates into sigmas with sigDe[]. It is exactly the same Dchi^2 that
   separates LCDM from the best fit in the (w0,wa) plane: that is why the
   significance matches how many sigmas of contour must be drawn to swallow LCDM.
   --------------------------------------------------------------------- *)

(* --- Helper to see the progress and time each fit ---
   PrintTemporary writes a notice that disappears when done (it tells you which
   one is running now); AbsoluteTiming measures the seconds. With the fast
   model (NDSolve) each fit takes seconds, not minutes.
   HoldRest prevents 'expr' from being evaluated before entering, so that it can
   be timed. *)
SetAttributes[ajusta, HoldRest];
ajusta[nombre_, expr_] := Module[{t, r},
   PrintTemporary["   ... ejecutando: ", nombre];
   {t, r} = AbsoluteTiming[expr];
   Print["   [OK] ", nombre, "   (", Round[t, 0.1], " s)"];
   r];

(* It is essential to clear the fit variables before minimising.
   NMinimize needs them as free symbols. If another file that sets h was evaluated
   before in the same kernel (for example h=0.7 in the LCDM plots one), that
   h bound to 0.7 sneaks into the list of variables ({Om,h,wb} -> {Om,0.7,wb})
   and produces 'NMinimize::ivar: 0.7 is not a valid variable' and the cascade of
   'ReplaceAll::reps'. This ClearAll makes the fit immune to the kernel state. *)
ClearAll[Om, h, w0, wa, wb];

(* --- (1) BAO --- *)
fitBAOl = ajusta["BAO / LCDM",    Quiet@NMinimize[{chi2BAO[Om, -1, 0], 0.1 < Om < 0.6}, Om]]; (* Gives values w0=-1, wa=0 and minimises over Om,
which is bounded between 0.1 and 0.6; LambdaCDM was known to be around 0.31, so it is a good bound. So we find min chi^2BAO in LCDM. Remember 
that H_0*r_d was already minimised analytically in BAO and we put its optimal value into the BAO formula.*)
fitBAOw = ajusta["BAO / w0waCDM", Quiet@NMinimize[{chi2BAO[Om, w0, wa],
     0.1 < Om < 0.6 && -3 < w0 < 1 && -3 < wa < 2 && w0 + wa < 0}, {Om, w0, wa}]]; (*The same for w0waCDM in the BAO-only case.*)
DchiBAO = First[fitBAOl] - First[fitBAOw];   (* Dchi2_MAP (>=0).            *)

(* A note on convergence with free wb.
   With wb (scale about 0.022) mixed with w0,wa (about 1) the default NMinimize
   gets stuck (it may return a w0waCDM worse than LCDM, that is Dchi2<0, which is
   absurd). It is fixed by forcing differential evolution, which spreads the
   search over the whole box of each parameter. The options go inline in each
   NMinimize, not as a variable, because NMinimize has HoldAll and would not
   evaluate them.                                                             *)

(* --- (2) BAO+BBN  (wb is now a free parameter with its BBN prior) --- *)
(* The bound 0.019<wb<0.025 is only the search box (about +-5 sigma); what
   really constrains wb is the Gaussian prior inside chi2Full.                *)
fitBBNl = ajusta["BAO+BBN / LCDM",    Quiet@NMinimize[{chi2Full[Om, h, -1, 0, wb, "bbn"],
     0.1 < Om < 0.6 && 0.5 < h < 0.9 && 0.019 < wb < 0.025}, {Om, h, wb},
     Method -> {"DifferentialEvolution", "SearchPoints" -> 60}, MaxIterations -> 500]]; (*chi^2min for LCDM in BAO+BBN and it also stores
     which optimal parameters give that minimum. *)
fitBBNw = ajusta["BAO+BBN / w0waCDM", Quiet@NMinimize[{chi2Full[Om, h, w0, wa, wb, "bbn"],
     0.1 < Om < 0.6 && 0.5 < h < 0.9 && 0.019 < wb < 0.025 &&
     -3 < w0 < 1 && -3 < wa < 2 && w0 + wa < 0}, {Om, h, w0, wa, wb},
     Method -> {"DifferentialEvolution", "SearchPoints" -> 60}, MaxIterations -> 500]]; (*chi^2min for w0waCDM in BAO+BBN and it also stores
     which optimal parameters give that minimum. *)
DchiBBN = First[fitBBNl] - First[fitBBNw]; (*Deltachi_C^2 BAO+BBN.*)

(* --- (3) BAO+CMB  (wb free; its prior and cross terms are in the 3x3 C). --- *)
fitCMBl = ajusta["BAO+CMB / LCDM",    Quiet@NMinimize[{chi2Full[Om, h, -1, 0, wb, "cmb"],
     0.1 < Om < 0.6 && 0.5 < h < 0.9 && 0.021 < wb < 0.024}, {Om, h, wb},
     Method -> {"DifferentialEvolution", "SearchPoints" -> 60}, MaxIterations -> 500]]; (*chi^2min for LCDM in BAO+CMB and it also stores
     which optimal parameters give that minimum.*)
fitCMBw = ajusta["BAO+CMB / w0waCDM", Quiet@NMinimize[{chi2Full[Om, h, w0, wa, wb, "cmb"],
     0.1 < Om < 0.6 && 0.5 < h < 0.9 && 0.021 < wb < 0.024 &&
     -3 < w0 < 1 && -3 < wa < 2 && w0 + wa < 0}, {Om, h, w0, wa, wb},
     Method -> {"DifferentialEvolution", "SearchPoints" -> 60}, MaxIterations -> 500]]; (*chi^2min for w0waCDM in BAO+CMB.*)
DchiCMB = First[fitCMBl] - First[fitCMBw]; (*Deltachi_C^2 BAO+CMB.*)

(* --- (4) BAO+CMB+SNe  (same scheme as BAO+CMB, with the SNe term). --- *)
fitCSl = ajusta["BAO+CMB+SNe / LCDM",    Quiet@NMinimize[{chi2CS[Om, h, -1, 0, wb],
     0.1 < Om < 0.6 && 0.5 < h < 0.9 && 0.021 < wb < 0.024}, {Om, h, wb},
     Method -> {"DifferentialEvolution", "SearchPoints" -> 60}, MaxIterations -> 500]]; (*chi^2min LCDM in BAO+CMB+SNe*)
fitCSw = ajusta["BAO+CMB+SNe / w0waCDM", Quiet@NMinimize[{chi2CS[Om, h, w0, wa, wb],
     0.1 < Om < 0.6 && 0.5 < h < 0.9 && 0.021 < wb < 0.024 &&
     -3 < w0 < 1 && -3 < wa < 2 && w0 + wa < 0}, {Om, h, w0, wa, wb},
     Method -> {"DifferentialEvolution", "SearchPoints" -> 60}, MaxIterations -> 500]]; (*chi^2min w0waCDM in BAO+CMB+SNe.*)
DchiCS = First[fitCSl] - First[fitCSw]; (*Deltachi_C^2 BAO+CMB+SNe.*)

Print["------------------------------------------------------------"];
Print["CASO      Dchi2_MAP     significancia (w0waCDM vs LCDM)"];
Print["BAO       ", NumberForm[DchiBAO, {4,2}], "        ", NumberForm[sigDe[DchiBAO], {3,2}], " sigma"];
Print["BAO+BBN   ", NumberForm[DchiBBN, {4,2}], "        ", NumberForm[sigDe[DchiBBN], {3,2}], " sigma"];
Print["BAO+CMB   ", NumberForm[DchiCMB, {4,2}], "        ", NumberForm[sigDe[DchiCMB], {3,2}], " sigma"];
Print["BAO+CMB+SNe ", NumberForm[DchiCS, {4,2}], "      ", NumberForm[sigDe[DchiCS], {3,2}], " sigma"]; 
Print["------------------------------------------------------------"];


(* =====================================================================
   Block 3b : common format of the plots (frame, grid, etc.)
   ===================================================================== *)
SetOptions[{Plot, ListContourPlot, ContourPlot, ListPlot, Graphics},
   Frame -> True, Axes -> False, GridLines -> Automatic,
   GridLinesStyle -> Directive[GrayLevel[0.85], AbsoluteThickness[0.4]],
   FrameStyle  -> Directive[Black, AbsoluteThickness[1]],
   LabelStyle  -> Directive[Black, 12], Background -> White,
   ImageSize   -> 440,
   ImagePadding -> {{Automatic, 24}, {Automatic, Automatic}}];


(* =====================================================================
   Block 4 : confidence regions (w0,wa) with LCDM marked
   =====================================================================
   How I get the confidence regions.
   I want the region in the (w0,wa) plane. The other parameters (Om, and h in
   BAO+CMB) are nuisances: at each point (w0,wa) of a grid I eliminate them by
   profiling, that is, I minimise chi^2 over them keeping (w0,wa) fixed:
        chi2_perf(w0,wa) = min_{Om,h}  chi^2(Om,h,w0,wa).
   Then I subtract the global minimum (the minimum over all the parameters,
   already computed before): Dchi^2(w0,wa) = chi2_perf - chi2_min.
   The X% confidence region is { (w0,wa) : Dchi^2 <= threshold(X) }.
   For two parameters the thresholds (Wilks, 2 d.o.f.) are:
        Dchi^2 = 2.30  -> 68.3% (about 1 sigma)
        Dchi^2 = 6.18  -> 95.4% (about 2 sigma)
        Dchi^2 = 11.83 -> 99.7% (about 3 sigma)
   (they come from inverting CL=1-e^{-Dchi^2/2}). I draw the contour lines of
   Dchi^2 at those values: they are the confidence contours.
   --------------------------------------------------------------------- *)

(* Profilers: they return chi^2 minimised over the nuisances, at every point
   (w0,wa) of the grid that we plot afterwards. This is not like the fits above,
   because there we simply minimised chi^2 with respect to all the parameters.
   Why 60 points and not 16: before, I profiled Om with only 16 values.
   The BAO degeneracy valley is long and almost flat (Dchi^2 between 0 and 2 over
   a huge stretch), so 16 points plus a coarse grid made the floor of the
   valley cross the 2.30 threshold erratically and the 68% region came out
   broken into islands. They were numerical islands, not physics. With a fine
   scan of Om (60 points) and more grid resolution, the profile is smooth and
   the region comes out connected.                                            *)
OmMalla = Subdivide[0.10, 0.60, 60];   (* Fine scan to profile Om.          *)
(* BAO: profiling Om is cheap (A is analytic), so Min over a fine scan. *)
profBAO[w0_, wa_] := If[w0 + wa >= 0, 60., Min[chi2BAO[#, w0, wa] & /@ OmMalla]]; (*To the points (w0,wa) s.t. w0+wa>=0 we assign chi^2=60,
which is a very large penalty (so that they do not fall in any confidence region). If the sum is negative, we evaluate the bao chi^2 over
the 60 Omega_m points we have generated and take the minimum of them.*)
(* BAO+CMB: profile (Om,h,wb) with local FindMinimum. A moderate AccuracyGoal
   means fewer iterations and therefore more speed, without losing sharpness. *)
profCMB[w0_, wa_] := If[w0 + wa >= 0, 60.,
   Quiet@First@FindMinimum[{chi2Full[Om, h, w0, wa, wb, "cmb"],
        0.1 < Om < 0.6 && 0.5 < h < 0.9 && 0.021 < wb < 0.024},
      {{Om, 0.35}, {h, 0.64}, {wb, 0.02236}}, AccuracyGoal -> 5, PrecisionGoal -> 5]]; (* As before, we penalise if the sum is positive.
      Then we look for a local minimum around a seed that we give to the parameters, because we know it is close to there, which makes
      the method much faster since here there are 3 free parameters, instead of the single one in BAO.*)

(* BAO+BBN: profile (Om,h,wb) with FindMinimum. It is fast because the ODE only
   goes up to z=2.5 (there is no z_star). Its contour comes out almost the same
   as the BAO one. *)
profBBN[w0_, wa_] := If[w0 + wa >= 0, 60.,
   Quiet@First@FindMinimum[{chi2Full[Om, h, w0, wa, wb, "bbn"],
        0.1 < Om < 0.6 && 0.5 < h < 0.9 && 0.019 < wb < 0.025},
      {{Om, 0.30}, {h, 0.68}, {wb, 0.02218}}, AccuracyGoal -> 5, PrecisionGoal -> 5]];(* As before, we penalise if the sum is positive.
      Then we look for a local minimum around a seed that we give to the parameters, because we know it is close to there, which makes
      the method much faster since here there are 3 free parameters, instead of the single one in BAO.*)

(* How to see the progress.
   I wrap each grid in  Monitor[ table , liveIndicator ]: while it runs,
   Mathematica shows a bar with the current node, so you see how far it has got
   and whether it stops at some point. AbsoluteTiming gives the total time. A
   global counter 'prog' counts the nodes. (To stop it: menu Evaluation >
   Abort Evaluation, or Alt+. on Windows / Cmd+. on Mac.)                     *)

(* ---- BAO grid (wide window: LCDM falls inside the 95%). ---- *)
w0aB = -1.6; w0bB = 0.4; waaB = -3.0; wabB = 1.5; nB = 34; (*For the plot we study w0, wa, -1.6<=w0<=0.4 and -3<=wa<=1.5 and on each axis, 34 interv. *)
w0sB = Subdivide[w0aB, w0bB, nB]; wasB = Subdivide[waaB, wabB, nB];
prog = 0; total = (nB + 1)^2;
{tBAO, gridBAO} = AbsoluteTiming@Monitor[
   Table[prog++; profBAO[w0, wa], {wa, wasB}, {w0, w0sB}], (*Table gives us a 35x35 matrix with the chi^2 minima at each point (w0,wa),
   and it is stored in gridBAO.*)
   Row[{"malla BAO:  nodo ", prog, " / ", total, "   ",
        ProgressIndicator[prog, {0, total}]}]];(*The progress indicator tells us in real time which iteration it is at out of the 1225 (35^2)
        it has to compute. AbsoluteTiming gives us, once all have been computed, the total time taken (it is stored in tBAO).*)
gridBAO = gridBAO - Min[gridBAO];(*Computes Deltachi_c^2 at each point, subtracting at each (w0,wa) gridBAO-min(all the gridBAO),
which is taken as the minimum of chi^2 over all the parameters, including w0,wa. Note that this minimum is over a discrete grid;
we could use fitBAOw, etc. which we defined before, but since this is only for plotting and the difference is so small, it is not noticeable.
It is fine to leave it like this. For the calculations of numerical values of sigmas and so on it has been done exactly.*)
Print["malla BAO lista en ", Round[tBAO, 0.1], " s"];(*Prints: "I have finished scanning the 1225 possible universes in X seconds". *)

(* ---- BAO+BBN grid (same wide window as BAO; it comes out almost the same). ---- *)
prog = 0; total = (nB + 1)^2;
{tBBN, gridBBN} = AbsoluteTiming@Monitor[
   Table[prog++; profBBN[w0, wa], {wa, wasB}, {w0, w0sB}],
   Row[{"malla BAO+BBN:  nodo ", prog, " / ", total, "   ",
        ProgressIndicator[prog, {0, total}]}]];
gridBBN = gridBBN - Min[gridBBN];
Print["malla BAO+BBN lista en ", Round[tBBN, 0.1], " s"]; (*Same as before.*)

(* ---- BAO+CMB grid (window centred between LCDM and the best fit) ----
   Here LCDM falls outside the 95%, so the window includes (-1,0) and the optimum
   (about -0.42,-1.73) and, to be able to see both and the 3 sigma contour, we
   recentre the limits of w0 and wa. It is the most expensive part, because it
   profiles Om,h with FindMinimum at each node. *)
w0aC = -1.25; w0bC = 0.15; waaC = -3.4; wabC = 0.9; nC = 34;
w0sC = Subdivide[w0aC, w0bC, nC]; wasC = Subdivide[waaC, wabC, nC];
prog = 0; total = (nC + 1)^2;
{tCMB, gridCMB} = AbsoluteTiming@Monitor[
   Table[prog++; profCMB[w0, wa], {wa, wasC}, {w0, w0sC}],
   Row[{"malla BAO+CMB:  nodo ", prog, " / ", total, "   ",
        ProgressIndicator[prog, {0, total}]}]];
gridCMB = gridCMB - Min[gridCMB];
Print["malla BAO+CMB lista en ", Round[tCMB, 0.1], " s"]; (*All the same.*)

(* ---------------------------------------------------------------------
   From here on everything is for drawing the plots. Contour drawer. It receives:
     grid      : the Dchi^2 grid,
     ventana   : {{w0min,w0max},{wamin,wamax}},
     niveles   : list of Dchi^2 thresholds to draw (2 or 3),
     sombras   : colours of the bands (one more than niveles),
     dLCDM     : exact Dchi^2 of LCDM (from the fits), for its label,
     offLab    : offset of the "LCDM" text relative to the point.
   It marks LCDM at (-1,0) with a red dot and a box with its exact
   confidence value (Dchi^2, sigmas and % CL).
   --------------------------------------------------------------------- *)
(* Translates a Dchi^2 threshold into its confidence label (for the legend). *)
etiquetaNivel[n_] := Which[
   Abs[n - 2.30]  < 0.1, "1\[Sigma] (68%)",
   Abs[n - 6.18]  < 0.1, "2\[Sigma] (95%)",
   Abs[n - 11.83] < 0.2, "3\[Sigma] (99.7%)",
   True, "nivel"];

(* --- Manual legend (fail-safe) ---
   SwatchLegend/LineLegend inside an Inset broke when drawn.
   Here I make the legend with text: a coloured square \[FilledSquare] plus its
   label, stacked in a Framed. Since it is text (like the LCDM box,
   which did come out fine), it is always drawn correctly.                    *)
swatch[col_] := Style["\[FilledSquare]", col, 13];
(* Legend glyphs are text only. A Graphics object placed inside the
   Framed/Inset legend rendered as a broken "double-click to edit" box (the very
   problem this code already noted for SwatchLegend/LineLegend), which made the
   star come out wrong and left a stray orange square. So the line and star
   samples are plain text characters, which always draw correctly.
   Coloured line-style text sample, matching the user's own notation:
   "solid" -> long dash (\[LongDash]), "dashed" -> "- -", "dotted" -> "\[CenterDot]\[CenterDot]\[CenterDot]". *)
lineTxt[col_, tipo_] := Style[
   tipo /. {"solid" -> "\[LongDash]", "dashed" -> "- -",
            "dotted" -> "\[CenterDot]\[CenterDot]\[CenterDot]"}, col, Bold, 12];
(* Gold star for the "(w0,wa) optimo" legend row (dark gold, so it is visible on
   the white legend background, unlike the bright yellow used on the dark fill). *)
estrellaTxt = Style["\[FivePointedStar]", RGBColor[0.85, 0.6, 0], 13];
(* Red dot for the LCDM legend row (same red as the LCDM point on the plot). *)
puntoLCDM = Style["\[FilledCircle]", Red, 12];
(* Each legend row may start with a colour (and then it gets a band square, as
   before) or with an already built glyph (line or symbol); glifoLeg picks which. *)
glifoLeg[x_] := If[MatchQ[x, _RGBColor | _GrayLevel | _Hue | _CMYKColor], swatch[x], x];
cajaLeg[titulo_, filas_] := Framed[
   Column[Prepend[(Row[{glifoLeg[#1], " ", #2}] &) @@@ filas, Style[titulo, 9, Bold]],
      Spacings -> 0.3, Alignment -> Left],
   Background -> White, FrameStyle -> Gray, RoundingRadius -> 3];
(* Solid colours of each band (for the small squares of the legend). *)
colBanda = {RGBColor[0.28, 0.44, 0.56],
            RGBColor[0.615, 0.7635, 0.857],
            RGBColor[0.829, 0.901, 0.955]};

(* Contour line styles, one per level, so that the levels are told apart also in
   black and white or by colour-blind readers, besides the fill shade: 68% solid,
   95% dashed, 99.7% dotted. Take[..,n] according to how many levels the figure
   draws (2 or 3). *)
estilosContorno[n_] := Take[{
    Directive[Black, AbsoluteThickness[1.1], Dashing[None]],           (* 68%   solid.  *)
    Directive[Black, AbsoluteThickness[1.1], AbsoluteDashing[{7, 4}]], (* 95%   dashed. *)
    Directive[Black, AbsoluteThickness[1.1], AbsoluteDashing[{1, 4}]]  (* 99.7% dotted. *)
   }, n];

(* Marker of the optimum (w0,wa) in the 4 individual figures: a yellow star
   with a black edge; it shows on the dark fill at the centre of the region and,
   by its shape, also in black and white (unlike the round LCDM dot). *)
estrellaOpt[{x_, y_}] := {Text[Style["\[FivePointedStar]", 18, Black], {x, y}],
                          Text[Style["\[FivePointedStar]", 15, RGBColor[1, 0.85, 0]], {x, y}]};

(* Generic optimum marker for the overlay plot: a glyph of a different shape
   per model (star/triangle/diamond), with a white halo so that it shows over the
   lines, in the colour of the model. *)
markOpt[{x_, y_}, glifo_, col_, tam_ : 14] := {
   Text[Style[glifo, tam + 3, White], {x, y}],
   Text[Style[glifo, tam, col], {x, y}]};

(* Axis box that encloses the painted region (Dchi^2 <= nivelMax) with a
   margin, clamped to the grid window (so that it never shows an area without
   data) and forcing the 'extra' points (LCDM and the optimum) in. This way each
   figure fits its axes to what is painted, with no dead space. Grid index: row
   -> wa, column -> w0. *)
cajaRegion[gf_, {{w0a_, w0b_}, {waa_, wab_}}, nivelMax_, extra_List, margen_ : 0.12] :=
  Module[{nr, nc, was, w0s, pos, xs, ys, mw0, mwa},
   {nr, nc} = Dimensions[gf];
   was = Subdivide[waa, wab, nr - 1];   (* Row i -> wa.    *)
   w0s = Subdivide[w0a, w0b, nc - 1];   (* Column j -> w0. *)
   pos = Position[gf, v_ /; v <= nivelMax, {2}];
   xs = Join[w0s[[pos[[All, 2]]]], extra[[All, 1]]];
   ys = Join[was[[pos[[All, 1]]]], extra[[All, 2]]];
   (* Tighter, cell-based margin so the region fills the frame. The old
      12% margin left a big empty band for the elongated BAO/BBN/CMB valleys (they
      looked "not adapted"). Now: at least ~1.5 grid cells (covers the linear-
      interpolation overshoot, so no painted region is ever cut) or 4% of the
      region span, whichever is larger. *)
   mw0 = Max[0.04 (Max[xs] - Min[xs]), 1.5 (w0b - w0a)/(nc - 1)];
   mwa = Max[0.04 (Max[ys] - Min[ys]), 1.5 (wab - waa)/(nr - 1)];
   {{Max[Min[xs] - mw0, w0a], Min[Max[xs] + mw0, w0b]},
    {Max[Min[ys] - mwa, waa], Min[Max[ys] + mwa, wab]}}];

(* (w0,wa) of the minimum of a grid (best fit on the grid): it always falls
   inside the window of that grid. Index: row -> wa, column -> w0. *)
optMalla[grid_, {{w0a_, w0b_}, {waa_, wab_}}] := Module[{im, d},
   d = Dimensions[grid]; im = First@Position[grid, Min[grid], {2}];
   {Subdivide[w0a, w0b, Last[d] - 1][[im[[2]]]],
    Subdivide[waa, wab, First[d] - 1][[im[[1]]]]}];

(* The optimum (w0,wa) is taken from the minimum of the plotted grid, so that
   it always falls inside the drawn window and right at the centre of the region;
   the BAO valley is so flat that the optimum of the global fit can land outside
   the window. It is marked with a star and added to the legend; the 'titulo'
   (name of the dataset) becomes the header of the legend, above it. *)
dibujaContorno[grid_, ventana_, niveles_, sombras_, titulo_, dLCDM_, offLab_, posLeg_] :=
  Module[{gf, w0a, w0b, waa, wab, itp, caja, optGrid, filasLeg},
   {{w0a, w0b}, {waa, wab}} = ventana;
   gf = GaussianFilter[grid, 1.5]; gf = gf - Min[gf];   (* Light smoothing. *)
   (* I draw with ContourPlot of a linear interpolation of the grid, not with
      ListContourPlot. The reason is that the small white patch inside the dark
      region was a filling glitch of ListContourPlot. The linear interpolation
      does not overshoot (it does not invent islands) and ContourPlot samples
      finely, so the filling comes out clean. Grid index: [wa,w0]. *)
   itp = ListInterpolation[gf, {{waa, wab}, {w0a, w0b}}, InterpolationOrder -> 1];
   optGrid = optMalla[grid, ventana];   (* Optimum (w0,wa) = minimum of the grid. *)
   (* Adapted axes: box enclosing the painted region (outermost level) plus LCDM and the optimum. *)
   caja = cajaRegion[gf, ventana, Max[niveles], {{-1, 0}, optGrid}];
   (* Legend rows (all TEXT glyphs): one row per level showing its band square +
      its contour line style + how many sigmas it is, then the optimum star, and
      last the LCDM red dot. *)
   filasLeg = Join[
      Table[{Row[{swatch[colBanda[[k]]], " ",
                  lineTxt[Black, {"solid", "dashed", "dotted"}[[k]]]}],
             etiquetaNivel[niveles[[k]]]}, {k, Length[niveles]}],
      {{estrellaTxt, Row[{"(", Subscript["w", "0"], ",", Subscript["w", "a"], ") \[OAcute]ptimo"}]}},
      {{puntoLCDM, "\[CapitalLambda]CDM"}}];
   ContourPlot[itp[wa, w0], {w0, w0a, w0b}, {wa, waa, wab},
     Contours      -> niveles,
     ContourShading-> sombras,
     ContourStyle  -> estilosContorno[Length[niveles]],  (* 68% solid, 95% dashed, 99.7% dotted (black and white / colour-blind). *)
     PlotPoints    -> 60, MaxRecursion -> 2,
     FrameLabel    -> {Style[Subscript["w", "0"], 13], Style[Subscript["w", "a"], 13]},
     PlotRange     -> caja,           (* Axes fitted to the painted region. *)
     PlotRangePadding -> None,
     Epilog -> {
        (* LCDM point (it is a point: the constant w=-1 model; the size of the
           dot is only cosmetic, it means nothing). *)
        Red, PointSize[0.028], Point[{-1, 0}],
        Text[Style["\[CapitalLambda]CDM", Red, Bold, 12], {-1, 0}, offLab],
        estrellaOpt[optGrid],   (* Optimum (w0,wa) of the fit, marked with a star. *)
        (* Manual legend: the name of the dataset is the header (above), then the
           confidence levels and, last, the optimum. *)
        Inset[cajaLeg[titulo, filasLeg], Scaled[posLeg]],   (* header = dataset; then LCDM, the sigma levels with their line style, and the optimum *)
        (* Box with the exact confidence value of LCDM. *)
        Inset[Framed[
           Style[Row[{"\[CapitalLambda]CDM:  \[CapitalDelta]",
                Superscript["\[Chi]", "2"], " = ", NumberForm[dLCDM, {3,1}],
                "  \[Rule]  ", NumberForm[sigDe[dLCDM], {2,1}], "\[Sigma]   (",
                NumberForm[100 clDe[dLCDM], {3,1}], "% CL)"}], 10, Black],
           Background -> White, FrameStyle -> Gray, RoundingRadius -> 3],
           Scaled[{0.5, 0.07}]]}]];

(* Band palettes (inside -> outside): dark blue, medium, light, white. *)
sombras2 = {RGBColor[0.28, 0.44, 0.56],       
            RGBColor[0.615, 0.7635, 0.857],    
            White};

sombras3 = {RGBColor[0.28, 0.44, 0.56],
            RGBColor[0.615, 0.7635, 0.857],
            RGBColor[0.829, 0.901, 0.955],     
            White};

(* --- BAO: LCDM is inside the 95%, so 68% and 95% are enough. --- *)
figBAO = dibujaContorno[gridBAO, {{w0aB, w0bB}, {waaB, wabB}},
   {2.30, 6.18}, sombras2, "BAO", DchiBAO, {1.5, -1.1}, {0.80, 0.82}];

(* --- BAO+BBN: almost the same as BAO (LCDM inside the 95%); I include it for
       completeness, to see where LCDM falls. --- *)
figBBNcont = dibujaContorno[gridBBN, {{w0aB, w0bB}, {waaB, wabB}},
   {2.30, 6.18}, sombras2, "BAO+BBN", DchiBBN, {1.5, -1.1}, {0.80, 0.82}];

(* --- BAO+CMB: LCDM is outside the 95%, so I add the 3 sigma (99.7%) to show
       that it does fall within 3 sigma. --- *)
figCMB = dibujaContorno[gridCMB, {{w0aC, w0bC}, {waaC, wabC}},
   {2.30, 6.18, 11.83}, sombras3, "BAO+CMB (comprimido)", DchiCMB, {-1.6, 1.2}, {0.80, 0.84}];

figBAO
figBBNcont
figCMB

(* =====================================================================
   Block 4 (continued) : confidence region (w0,wa) of BAO+CMB+SNe
   ---------------------------------------------------------------------
   Same method as BAO+CMB (at each node (w0,wa) I profile Om,h,wb with
   FindMinimum), but with chi2CS. Adding the SNe narrows the region a lot
   and brings it closer to LCDM: the best fit lies near (w0,wa)~(-0.85,-0.50) and
   LCDM(-1,0) falls between 2 and 3 sigma, so I again use the three levels
   68/95/99.7 and a centred window smaller than the BAO+CMB one.               *)
profCS[w0_, wa_] := If[w0 + wa >= 0, 60.,
   Quiet@First@FindMinimum[{chi2CS[Om, h, w0, wa, wb],
        0.1 < Om < 0.6 && 0.5 < h < 0.9 && 0.021 < wb < 0.024},
      {{Om, 0.31}, {h, 0.675}, {wb, 0.02247}}, AccuracyGoal -> 5, PrecisionGoal -> 5]];
w0aS = -1.15; w0bS = -0.55; waaS = -1.8; wabS = 0.7; nS = 34; (* (w0,wa) window of the SNe case. *)
w0sS = Subdivide[w0aS, w0bS, nS]; wasS = Subdivide[waaS, wabS, nS];
prog = 0; total = (nS + 1)^2;
{tCS, gridCS} = AbsoluteTiming@Monitor[
   Table[prog++; profCS[w0, wa], {wa, wasS}, {w0, w0sS}],
   Row[{"malla BAO+CMB+SNe:  nodo ", prog, " / ", total, "   ",
        ProgressIndicator[prog, {0, total}]}]];
gridCS = gridCS - Min[gridCS];
Print["malla BAO+CMB+SNe lista en ", Round[tCS, 0.1], " s"];
figCS = dibujaContorno[gridCS, {{w0aS, w0bS}, {waaS, wabS}},
   {2.30, 6.18, 11.83}, sombras3, "BAO+CMB+SNe (comprimido)", DchiCS, {-1.5, 1.1}, {0.80, 0.84}];
figCS


(* =====================================================================
   Block 4b : H0 and the Hubble tension (the six DESI fits)
   =====================================================================
   The key idea is that H0 here is not measured, it is inferred, and what fixes
   its value is the model, not the dataset:
     - In LCDM (w=-1, rigid) the three datasets, BAO+BBN, BAO+CMB and
       BAO+CMB+SNe, give the same H0 (about 68.5): the BAO shape fixes Om and the
       ruler rd (calibrated by BBN or by CMB, which agree) fixes the scale, so
       there is a single H0. The three LCDM bars almost overlap.
     - In w0waCDM the (w0,wa) freedom opens the geometric degeneracy and H0
       drops to about 63 with BAO+BBN and BAO+CMB; but adding the SNe
       (BAO+CMB+SNe) H0 rises to about 67.6, because the SNe bring (w0,wa) closer
       to LCDM and reduce the degeneracy, which is why that bar sits next to the
       LCDM group.
   That is why I draw all six fits (3 datasets x 2 models) grouped by
   model (blues = LCDM, purples = w0waCDM), plus two external anchors: SH0ES
   (local, late) and Planck (LCDM, early). This shows that the jump from about
   68.5 to about 63 is due to changing the model, not to changing the data.

   The error bars (the sigma of each H0) are not set by hand: the code computes
   them by a one-dimensional profile (the same method as Block 4c, but on the H0
   axis). For each fixed H0 I minimise chi^2 over the rest of the parameters and
   find where Dchi^2 rises to 1 (the Wilks threshold for 1 d.o.f. = 1 sigma); the
   half-width between the two crossings is the sigma. This is done by
   profChiH0[] and errH0[], defined below. Only SH0ES (1.04) and Planck (0.54)
   have a fixed sigma: they are published values, not my fits.
   --------------------------------------------------------------------- *)
(* ---------------------------------------------------------------------
   One-dimensional profile of H0 (the error bars are computed here)
   ---------------------------------------------------------------------
   profChiH0[caso, modelo, hf, semilla] fixes H0 = 100*hf and returns the
   lowest possible chi^2 moving the rest of the parameters for that fixed h.
   'semilla' are the rules of the global best fit; I start FindMinimum from there
   so that it tracks the minimum stably as hf moves away, and does not get lost,
   as did happen to NMinimize with free wb. 'caso' is "bbn" or "cmb"; 'modelo'
   decides which parameters remain free:
     - "lcdm" : {Om, wb} are profiled           (w0=-1, wa=0 fixed).
     - "w0wa" : {Om, w0, wa, wb} are profiled   (dynamical dark energy).
   The bounds are only the local search box; the real prior on wb already
   lives inside chi2Full (Gaussian BBN / CMB inside the 3x3 C).               *)
profChiH0[caso_, "lcdm", hf_?NumericQ, semilla_] := Quiet@First@FindMinimum[
   {chi2Full[Om, hf, -1, 0, wb, caso], 0.1 < Om < 0.6 && 0.019 < wb < 0.025},
   {{Om, Om /. semilla}, {wb, wb /. semilla}}];
profChiH0[caso_, "w0wa", hf_?NumericQ, semilla_] := Quiet@First@FindMinimum[
   {chi2Full[Om, hf, w0, wa, wb, caso],
    0.1 < Om < 0.6 && -3 < w0 < 1 && -3 < wa < 2 && w0 + wa < 0 && 0.019 < wb < 0.025},
   {{Om, Om /. semilla}, {w0, w0 /. semilla}, {wa, wa /. semilla}, {wb, wb /. semilla}}];
(* The same two definitions of the one-dimensional profile of H0, but for
   BAO+CMB+SNe (they use chi2CS). Since they carry the literal "cs", Mathematica
   prefers them over the generic ones above (caso_) when calling
   profChiH0["cs", ...], so errH0 works the same without touching it.          *)
profChiH0["cs", "lcdm", hf_?NumericQ, semilla_] := Quiet@First@FindMinimum[
   {chi2CS[Om, hf, -1, 0, wb], 0.1 < Om < 0.6 && 0.021 < wb < 0.024},
   {{Om, Om /. semilla}, {wb, wb /. semilla}}];
profChiH0["cs", "w0wa", hf_?NumericQ, semilla_] := Quiet@First@FindMinimum[
   {chi2CS[Om, hf, w0, wa, wb],
    0.1 < Om < 0.6 && -3 < w0 < 1 && -3 < wa < 2 && w0 + wa < 0 && 0.021 < wb < 0.024},
   {{Om, Om /. semilla}, {w0, w0 /. semilla}, {wa, wa /. semilla}, {wb, wb /. semilla}}];

(* errH0 turns that profile into the 1 sigma bar:
     1) it scans h on a fine grid around the best fit,
     2) it computes Dchi^2(h) = profChiH0(h) - chi2min; we look for the h above
        and below where Deltachi^2=1 (according to Wilks this is the 1 sigma case
        for a single parameter, which here is h),
     3) it locates the two h where Dchi^2 crosses 1 (one below and one
        above the minimum), because we are on a discrete grid, and interpolates
        the crossing linearly (the profile is a smooth parabola near the minimum,
        so the straight line between two very close neighbouring points is enough,
        it is a first-order Taylor expansion),
     4) it returns the half-width in H0 = 100*(h_right - h_left)/2, which is the
        sigma.
   'ic' = index of the minimum of the profile; SelectFirst looks, on each side,
   for the first point that already exceeds Dchi^2=1; 'cruza[i,j]' interpolates
   the h where Dchi^2=1 between that point (j) and the previous one (i, where it
   was still below 1). *)
errH0[caso_, modelo_, chi2min_, hbest_?NumericQ, semilla_] := Module[
   {paso = 0.001, hs, dch, ic, cruza, kDer, kIzq},
   hs  = Range[hbest - 0.08, hbest + 0.08, paso];                  (* h grid.                 *)
   dch = (profChiH0[caso, modelo, #, semilla] - chi2min) & /@ hs;  (* Dchi^2(h) on the grid.  *)
   ic  = First@Ordering[dch, 1];                                   (* lowest point (minimum) of the chi^2 parabola (it takes the minimum grid 
   point, not the vertex of the parabola, although they are very similar, and the difference is negligible *)
   cruza[i_, j_] := hs[[i]] + (1 - dch[[i]]) (hs[[j]] - hs[[i]])/(dch[[j]] - dch[[i]]);(*I define the interpolation between any two values of h,
   neighbouring ones.*)
   kDer = SelectFirst[Range[ic, Length[hs]], dch[[#]] >= 1 &];     (* First crossing (right).  *)
   kIzq = SelectFirst[Range[ic, 1, -1],      dch[[#]] >= 1 &];     (* First crossing (left).   *)
   100 (cruza[kDer - 1, kDer] - cruza[kIzq + 1, kIzq])/2];         (*I take h of the right point - h of the left point = 2 sigma (h).
    I solve for sigma(h) and multiply by 100 so that it is sigma (H0) [km/s/Mpc]. *)

(* Central H0 of each fit (H0 = 100 h; I read the optimal h from the rules of each fit). *)
h0BbnL = 100 (h /. fitBBNl[[2]]);   (* BAO+BBN , LCDM.     *)
h0CmbL = 100 (h /. fitCMBl[[2]]);   (* BAO+CMB , LCDM.     *)
h0BbnW = 100 (h /. fitBBNw[[2]]);   (* BAO+BBN , w0waCDM.  *)
h0CmbW = 100 (h /. fitCMBw[[2]]);   (* BAO+CMB , w0waCDM.  *)
h0CsL  = 100 (h /. fitCSl[[2]]);    (* BAO+CMB+SNe , LCDM.      *)
h0CsW  = 100 (h /. fitCSw[[2]]);    (* BAO+CMB+SNe , w0waCDM.   *)

(* One-dimensional sigma of each H0, computed by profiling. Arguments of errH0:
   caso, modelo, chi2_min of the fit (First@fit), h of the fit, and the rules
   of the fit as the seed to start the profiling.                             *)
eBbnL = errH0["bbn", "lcdm", First[fitBBNl], h /. fitBBNl[[2]], fitBBNl[[2]]];
eCmbL = errH0["cmb", "lcdm", First[fitCMBl], h /. fitCMBl[[2]], fitCMBl[[2]]];
eBbnW = errH0["bbn", "w0wa", First[fitBBNw], h /. fitBBNw[[2]], fitBBNw[[2]]];
eCmbW = errH0["cmb", "w0wa", First[fitCMBw], h /. fitCMBw[[2]], fitCMBw[[2]]];
eCsL  = errH0["cs",  "lcdm", First[fitCSl], h /. fitCSl[[2]], fitCSl[[2]]]; 
eCsW  = errH0["cs",  "w0wa", First[fitCSw], h /. fitCSw[[2]], fitCSw[[2]]]; 

datosH0 = {   (* {y, label, H0, sigma, colour}  (y: 1 bottom .. 8 top). *)
   {8, "SH0ES 2022 (local)",                   73.04,  1.04,  RGBColor[0.75,0.25,0.20]},
   {7, "DESI BAO+BBN (\[CapitalLambda]CDM)",   h0BbnL, eBbnL, RGBColor[0.10,0.30,0.45]},
   {6, "DESI BAO+CMB (\[CapitalLambda]CDM)",   h0CmbL, eCmbL, RGBColor[0.20,0.45,0.65]},
   {5, "DESI BAO+CMB+SNe (\[CapitalLambda]CDM)", h0CsL, eCsL, RGBColor[0.35,0.62,0.82]}, 
   {4, "Planck 2018 (\[CapitalLambda]CDM)",    67.36,  0.54,  RGBColor[0.20,0.55,0.35]},
   {3, "DESI BAO+BBN (w0waCDM)",               h0BbnW, eBbnW, RGBColor[0.55,0.35,0.65]},
   {2, "DESI BAO+CMB (w0waCDM)",               h0CmbW, eCmbW, RGBColor[0.45,0.25,0.55]},
   {1, "DESI BAO+CMB+SNe (w0waCDM)",           h0CsW, eCsW, RGBColor[0.70,0.45,0.85]}}; 

figBBN = Graphics[{
   Table[With[{y = datosH0[[i,1]], v = datosH0[[i,3]], e = datosH0[[i,4]],
               col = datosH0[[i,5]]},
     {col, AbsoluteThickness[2.4],
      Line[{{v - e, y}, {v + e, y}}],                                  (* Error bar. *)
      Line[{{v - e, y - 0.13}, {v - e, y + 0.13}}],                    (* Left cap.  *)
      Line[{{v + e, y - 0.13}, {v + e, y + 0.13}}],                    (* Right cap. *)
      PointSize[0.024], Point[{v, y}],                                 (* Central value. *)
      Black, Text[Style[Row[{NumberForm[v, {4,2}], " \[PlusMinus] ", NumberForm[e, {3,2}]}], 10],
                  {v + e, y}, {-1.2, 0}]}],                            (* Value on the right. *)
     {i, Length[datosH0]}]},
   Frame -> True, Axes -> False,
   FrameLabel -> {Row[{Subscript["H", "0"], "   [km/s/Mpc]"}], ""},
   FrameTicks -> {{Table[{datosH0[[i,1]], Style[datosH0[[i,2]], 10]}, {i, Length[datosH0]}], None},
                  {Automatic, None}},
   PlotRange -> {{57, 81}, {0.5, 8.5}},
   GridLines -> {Range[58, 80, 2], None},
   GridLinesStyle -> Directive[GrayLevel[0.85], AbsoluteThickness[0.4]],
   AspectRatio -> 0.68, ImageSize -> 640, Background -> White];
figBBN


(* =====================================================================
   Block 4c : constraint on Omega_M (one-dimensional likelihood, LCDM)
   =====================================================================
   What it is for: to see how well each combination measures the matter
   density Om, and how the CMB tightens it. It is a one-dimensional profile: for
   each fixed Om, I minimise chi^2 over everything else (LCDM: w0=-1, wa=0 fixed;
   h is profiled; in BAO the scale A is analytic).
   Dchi^2(Om) = chi2_perf(Om) - chi2_min.

   I plot it as a likelihood  L(Om) = e^{-Dchi^2/2}  (normalised to peak 1),
   which reads as a "probability distribution" of Om (under a flat prior,
   proportional to the posterior). It is more intuitive than the Dchi^2 parabola:
   the narrower the bell, the more precise the measurement (here, the CMB). The
   basis of all this: the pure absolute likelihood function (Likelihood,
   $\mathcal{L}$) is defined using the total $\chi^2$,
   not the $\Delta\chi^2$. That is:$$\mathcal{L} \propto e^{-\chi^2/2}$$However, if you try to plot that directly,
   you will run into a huge computational and visual problem. If your best model has a $\chi^2_{min} = 15$, its absolute probability
   would be $e^{-7.5} \approx 0.00055$. You would have to draw a tiny bell squashed against the X axis, impossible to read.
   To fix it, in cosmology (and in statistics in general) we always normalise the likelihood so that the peak of the bell equals
   exactly 1. How is it normalised mathematically? By dividing the whole curve by the maximum possible value
   ($\mathcal{L}_{max}$):$$\mathcal{L}_{norm} = \frac{\mathcal{L}}{\mathcal{L}_{max}} = \frac{e^{-\chi^2/2}}{e^{-\chi^2_{min}/2}}$$
   By the properties of exponents, when dividing powers with the same base, the exponents are subtracted:
   $$\mathcal{L}_{norm} = e^{-(\chi^2 - \chi^2_{min})/2} = e^{-\Delta\chi^2/2}$$There you have it!
   Using the $\Delta\chi^2$ is mathematically equivalent to plotting the absolute likelihood normalised to 1. The chi_^2min for each model
   is reached at one point of the whole grid.

   Confidence for a single parameter (different from the 2D contours).
   With only one parameter, the X% interval is where Dchi^2 <= threshold, with the
   Wilks thresholds for 1 d.o.f.:  1 (68%), 4 (95%), 9 (99.7%). In likelihood
   these are the lines  L = e^{-1/2}=0.61 (68%)  and  e^{-4/2}=0.14 (95%).
   Note that the interval is wider the more sure one wants to be, so more values
   of Om are included; that is why the 95% cut is wider than the 68% one. Fixed
   parameters: w0,wa (that is, LCDM), wb, Neff, wnu.                          *)
profOmBAO[om_?NumericQ] := chi2BAO[om, -1, 0];
profOmBBN[om_?NumericQ] := Quiet@First@FindMinimum[
    {chi2Full[om, h, -1, 0, wb, "bbn"], 0.5 < h < 0.9 && 0.019 < wb < 0.025},
    {{h, 0.68}, {wb, 0.02218}}];
profOmCMB[om_?NumericQ] := Quiet@First@FindMinimum[
    {chi2Full[om, h, -1, 0, wb, "cmb"], 0.5 < h < 0.9 && 0.021 < wb < 0.024},
    {{h, 0.68}, {wb, 0.02236}}];
profOmCS[om_?NumericQ] := Quiet@First@FindMinimum[  (* Profile of Om with BAO+CMB+SNe (LCDM). *)
    {chi2CS[om, h, -1, 0, wb], 0.5 < h < 0.9 && 0.021 < wb < 0.024},
    {{h, 0.68}, {wb, 0.02236}}];

OmEje = Subdivide[0.26, 0.34, 120];
vero[ch_] := Exp[-(ch - Min[ch])/2];
LB = vero[profOmBAO /@ OmEje];(*gives us a vector of likelihood values L = e^{-Dchi^2/2}, normalised,
evaluated at the points of the Omega_m grid *)
LC = vero[profOmCMB /@ OmEje];(*gives us a vector of likelihood values L = e^{-Dchi^2/2}, normalised,
evaluated at the points of the Omega_m grid *)
LCS = vero[profOmCS /@ OmEje];(* One-dimensional likelihood of Om for BAO+CMB+SNe (LCDM). *)
(* Note: BAO+BBN gives exactly the same Om curve as BAO (BBN only fixes the
   scale rd/H0, it does not change how Om is measured). So I do not plot it
   apart, because one would hide the other: I put a single curve "BAO=BAO+BBN". *)

(* Om of the peak of each curve = the Om that maximises the one-dimensional
   likelihood (equivalent to the best-fit Om of each dataset). I take it from the
   already computed vector with Ordering (position of the maximum) so that the
   vertical line falls exactly on the plotted peak.                            *)
omPicoB  = OmEje[[First@Ordering[LB,  -1]]];   (* BAO (= BAO+BBN). *)
omPicoC  = OmEje[[First@Ordering[LC,  -1]]];   (* BAO+CMB.         *)
omPicoCS = OmEje[[First@Ordering[LCS, -1]]];   (* BAO+CMB+SNe.     *)

figOm = ListLinePlot[
   {Transpose[{OmEje, LB}], Transpose[{OmEje, LC}], Transpose[{OmEje, LCS}]}, (* Third curve: BAO+CMB+SNe. *)
   PlotStyle -> {Directive[RGBColor[0.20,0.50,0.80], AbsoluteThickness[2.2], Dashing[None]],
                 Directive[RGBColor[0.85,0.30,0.15], AbsoluteThickness[2.6], AbsoluteDashing[{9,5}]],
                 Directive[RGBColor[0.15,0.55,0.30], AbsoluteThickness[2.6], AbsoluteDashing[{2,4}]]}, (* Dash per curve on top of colour: BAO 
                 solid, BAO+CMB dashed, BAO+CMB+SNe dotted (black and white / colour-blind). *)
   Frame -> True, Axes -> False, GridLines -> {Automatic, None},
   GridLinesStyle -> Directive[GrayLevel[0.85], AbsoluteThickness[0.4]],
   FrameLabel -> {Subscript["\[CapitalOmega]", "M"], "verosimilitud  L  (pico = 1)"}, (* Omega_M axis with uppercase M. *)
   PlotRange -> {{0.26, 0.34}, {0, 1.08}},
   Epilog -> {
      GrayLevel[0.55], AbsoluteDashing[{4, 3}],
      Line[{{0.26, 0.607}, {0.34, 0.607}}], Line[{{0.26, 0.135}, {0.34, 0.135}}],
      Text[Style["68% (1\[Sigma])", 8, Gray], {0.338, 0.607}, {1, -0.6}],
      Text[Style["95% (2\[Sigma])", 8, Gray], {0.338, 0.135}, {1, -0.6}],
      (* Dash-dot vertical line from the peak of each curve down to the Om axis,
         with the value of Om (2 decimals) where it cuts. The colour is that of
         its curve, and the style is the same dash-dot as the era-transition
         vertical lines.                                                        *)
      {Directive[RGBColor[0.20,0.50,0.80], AbsoluteDashing[{6,3,1,3}], AbsoluteThickness[1.3]],
         Line[{{omPicoB, 0}, {omPicoB, 1}}]},
      {Directive[RGBColor[0.85,0.30,0.15], AbsoluteDashing[{6,3,1,3}], AbsoluteThickness[1.3]],
         Line[{{omPicoC, 0}, {omPicoC, 1}}]},
      {Directive[RGBColor[0.15,0.55,0.30], AbsoluteDashing[{6,3,1,3}], AbsoluteThickness[1.3]],
         Line[{{omPicoCS, 0}, {omPicoCS, 1}}]},
      Text[Framed[Style[NumberForm[omPicoB, {3,2}], 9, RGBColor[0.20,0.50,0.80]],
            Background -> White, FrameStyle -> None, FrameMargins -> 1], {omPicoB, 0.10}, {1.2, 0}],
      Text[Framed[Style[NumberForm[omPicoC, {3,2}], 9, RGBColor[0.85,0.30,0.15]],
            Background -> White, FrameStyle -> None, FrameMargins -> 1], {omPicoC, 0.10}, {-1.2, 0}],
      Text[Framed[Style[NumberForm[omPicoCS, {3,2}], 9, RGBColor[0.15,0.55,0.30]],
            Background -> White, FrameStyle -> None, FrameMargins -> 1], {omPicoCS, 0.24}, {-1.2, 0}],
      Inset[cajaLeg["perfil 1D (\[CapitalLambda]CDM)",
         {{lineTxt[RGBColor[0.20,0.50,0.80], "solid"],  "BAO  (= BAO+BBN)"},
          {lineTxt[RGBColor[0.85,0.30,0.15], "dashed"], "BAO+CMB"},
          {lineTxt[RGBColor[0.15,0.55,0.30], "dotted"], "BAO+CMB+SNe"}}], Scaled[{0.83, 0.82}]]}, (* Text line glyphs (no Graphics) so the 
          legend always renders *)
   ImageSize -> 500];
figOm


(* =====================================================================
   Block 4d : (w0,wa) overlay: BAO vs BAO+CMB vs BAO+CMB+SNe (key plot)
   =====================================================================
   Nothing new is computed here: the w0-wa plots of BAO, BAO+CMB and BAO+CMB+SNe
   are simply overlaid to see the difference.
   What it is for: it puts on the same axes the 68% and 95% regions of BAO
   (large: LCDM inside), of BAO+CMB (small: LCDM outside the 95%) and of
   BAO+CMB+SNe (even smaller: LCDM between 2 and 3 sigma). One sees at a glance
   how the CMB narrows the region and shifts it to the dynamical dark energy
   quadrant (w0>-1, wa<0), and how the SNe tighten it even more, bringing it
   back towards LCDM. There lies the tension with LCDM and the origin of the
   DESI result. It reuses gridBAO, gridCMB and gridCS, so it recomputes nothing.
   Confidence: 68% (2.30) and 95% (6.18), 2 d.o.f. Dashed = 68%; solid = 95%. *)
(* Three very different colours: BAO blue, BAO+CMB red, BAO+CMB+SNe green;
   within each colour the level is the line style (68% dashed, 95% solid), as in
   the compact legend the user prefers. The three optimum markers (star/triangle/
   diamond, one shape per model) give the extra black-and-white / colour-blind cue.
   I smooth the same way as the contours. *)
soloLineas[grid_, ventana_, col_] := Module[{gf, w0a, w0b, waa, wab, itp},
   {{w0a, w0b}, {waa, wab}} = ventana;
   gf = GaussianFilter[grid, 2]; gf = gf - Min[gf];
   itp = ListInterpolation[gf, {{waa, wab}, {w0a, w0b}}, InterpolationOrder -> 1];
   ContourPlot[itp[wa, w0], {w0, w0a, w0b}, {wa, waa, wab},
     Contours -> {2.30, 6.18}, ContourShading -> None, PlotPoints -> 60, MaxRecursion -> 2,
     ContourStyle -> {Directive[col, AbsoluteThickness[1.3], AbsoluteDashing[{9, 7}]], (* 68% dashed *)
                      Directive[col, AbsoluteThickness[2.8]]}]];                       (* 95% solid *)

(* Each model with its own dash (solid/dashed/dash-dot) on top of the
   colour, and its optimum (w0,wa) marked with a different shape (star / triangle
   / diamond). This way the three models and their optima are told apart also in
   black and white or by colour-blind readers. *)
figSolape = Show[
   soloLineas[gridBAO, {{w0aB, w0bB}, {waaB, wabB}}, RGBColor[0.15,0.45,0.80]],(* BAO         *)
   soloLineas[gridCMB, {{w0aC, w0bC}, {waaC, wabC}}, RGBColor[0.82,0.15,0.15]],(* BAO+CMB     *)
   soloLineas[gridCS, {{w0aS, w0bS}, {waaS, wabS}}, RGBColor[0.15,0.55,0.30]], (* BAO+CMB+SNe *)
   Frame -> True, Axes -> False, GridLines -> Automatic,
   GridLinesStyle -> Directive[GrayLevel[0.85], AbsoluteThickness[0.4]],
   FrameLabel -> {Style[Subscript["w", "0"], 13], Style[Subscript["w", "a"], 13]},
   PlotRange -> {{-1.6, 0.3}, {-3, 1}},   (* X axis unchanged; y axis from -3 to 1. *)
   Epilog -> {Red, PointSize[0.012], Point[{-1, 0}],
      Text[Style["\[CapitalLambda]CDM", Red, Bold, 12], {-1, 0}, {1.4, -1}],
      (* Optimum (w0,wa) of each model (minimum of each grid): a different
         shape per model, in the model colour, so they are told apart in B&W too. *)
      markOpt[optMalla[gridBAO, {{w0aB,w0bB},{waaB,wabB}}], "\[FivePointedStar]",  RGBColor[0.15,0.45,0.80]],
      markOpt[optMalla[gridCMB, {{w0aC,w0bC},{waaC,wabC}}], "\[FilledUpTriangle]", RGBColor[0.82,0.15,0.15]],
      markOpt[optMalla[gridCS,  {{w0aS,w0bS},{waaS,wabS}}], "\[FilledDiamond]",    RGBColor[0.15,0.55,0.30]],
      (* Compact legend (text only, so it always renders and fits): first
         line = what the dashed and solid lines mean; then one row per model
         (colour square + name); then the three optimum shapes lined up under the
         model colours + "(w0,wa) optimos" on the right, and last the LCDM red dot.
         Anchored by its top-right corner so it never spills out of the frame. *)
      Inset[cajaLeg["68% (- -)  y  95% (\[LongDash])",
         {{RGBColor[0.15,0.45,0.80], "BAO"},
          {RGBColor[0.82,0.15,0.15], "BAO+CMB"},
          {RGBColor[0.15,0.55,0.30], "BAO+CMB+SNe"},
          {Row[{Style["\[FivePointedStar]", RGBColor[0.15,0.45,0.80], 12], " ",
                Style["\[FilledUpTriangle]", RGBColor[0.82,0.15,0.15], 12], " ",
                Style["\[FilledDiamond]", RGBColor[0.15,0.55,0.30], 12]}],
           Row[{"(", Subscript["w","0"], ",", Subscript["w","a"], ") \[OAcute]ptimos"}]},
          {puntoLCDM, "\[CapitalLambda]CDM"}}],
         Scaled[{0.98, 0.97}], Scaled[{1, 1}]]},
   ImageSize -> 470];
figSolape


(* =====================================================================
   Block 5 : comparison tables with DESI (framed)
   ===================================================================== *)
estiloTabla[data_, cab_] := Grid[Prepend[data, Style[#, Bold] & /@ cab],
   Frame -> All, Alignment -> Left, Spacings -> {2, 0.8},
   Background -> {None, {1 -> RGBColor[0.80, 0.87, 0.92]}},   (* Blue header. *)
   BaseStyle -> 11];

(* Table A: parameters of each fit against DESI DR2, with their 1D sigma ---
   Each fit determines several parameters; I compare them with DESI putting one
   row per parameter. I split it into two tables by model:
     Table A1 = the LCDM cases      (w0=-1, wa=0 are fixed, so they are not compared).
     Table A2 = the w0waCDM cases   (BAO+CMB and BAO+CMB+SNe).
   The first column (the case) is extended with SpanFromAbove (it merges cells;
   the Grid uses Frame->All). Free parameters:
     BAO (LCDM)         -> Om, h*rd   (BAO alone only measures h*rd, not H0 and rd separately).
     BAO+BBN (LCDM)     -> Om, H0, wb.
     BAO+CMB (LCDM)     -> Om, H0, wb.
     BAO+CMB (w0wa)     -> Om, H0, w0, wa, wb.
     BAO+CMB+SNe (LCDM) -> Om, H0, wb.
     BAO+CMB+SNe (w0wa) -> Om, H0, w0, wa, wb.

   Why in w0waCDM I show BAO+CMB and BAO+CMB+SNe, and not BAO or BAO+BBN: in
   those two the model is degenerate. Without the CMB anchor at z~1100, w0 and wa
   do not separate (geometric degeneracy) and the best fit runs to the edge of the
   box (w0~-0.18, wa~-2.7). Moreover, there w0waCDM barely improves on
   LCDM (only about 1.7-1.9 sigma, Table B), which is not enough tension with
   LCDM. That is why DESI only publishes w0waCDM combined with CMB (or SNe), and
   so do I.

   Sigma: the "mi ajuste" (my fit) column is now value +- sigma. The 1D
   sigma of each parameter is computed by the code (function 'sigCurv',
   below): it profiles the parameter (fixes it, minimises chi^2 over the rest)
   and measures where Dchi^2 rises to 1 (the Wilks threshold for 1 d.o.f.). Note
   that for BAO+CMB in w0waCDM the sigma does not come out 0: the best fit is
   a point, but its uncertainty (how far it can move with Dchi^2<=1) is finite:
   it gives w0=-0.45+-0.25, wa=-1.64+-0.73, and so on, similar to the bars DESI
   publishes. (w0 and wa are very correlated; this sigma is the marginal one of
   each; their joint (w0,wa) ellipse is drawn in Block 4d.)

   DESI DR2 values (2503.14738), only the central reference value:
     eq.17 -> BAO (LCDM):        Om=0.2975 , h*rd=101.54
     eq.19 -> BAO+BBN (LCDM):    Om=0.2977 , H0=68.51 ; wb=0.02218 BBN prior
     Table 3 -> BAO+CMB (LCDM):  Om=0.301 , H0=68.4 ; wb=0.02236 (Planck)
     eq.24 -> BAO+CMB (w0waCDM): w0=-0.42, wa=-1.75, Om=0.353, H0=63.6, wb=0.02236
     DESI+CMB+Pantheon+ (w0waCDM): w0=-0.838, wa=-0.62, Om=0.3114, H0=67.51, wb=0.02236
   Note that my "CMB" is the compressed version (three numbers); it reproduces the
   full DESI+CMB approximately, which is why my values fall close but are not
   identical. *)

(* ---- Machinery for the 1D sigma of each parameter -----------------
   chi2BAOA: BAO chi^2 with the scale A explicit. In the rest of BAO, A is
   marginalised analytically; here I need it as a real variable
   to be able to measure how much h*rd = c/(100 A) moves.                      *)
chi2BAOA[Om_?NumericQ, A_?NumericQ] := Module[{orr = wr/0.49, Iz, g, res},
   Iz  = Ifun[Om, -1, 0, orr, 2.5];
   g   = gvecOf[Iz, Om, -1, 0, orr];
   res = dvec - A g;
   res . Cinv . res];

(* LCDM profiles (except BAO alone): I fix one parameter at 'v' and minimise over
   the other two (w0=-1, wa=0 fixed). 's' = best-fit rules, to start nearby.   *)
profL[caso_, "Om", v_, s_] := First@Quiet@FindMinimum[
   {chi2Full[v, h, -1, 0, wb, caso], 0.5 < h < 0.9 && 0.019 < wb < 0.025},
   {{h, h /. s}, {wb, wb /. s}}];
profL[caso_, "h", v_, s_] := First@Quiet@FindMinimum[
   {chi2Full[Om, v, -1, 0, wb, caso], 0.1 < Om < 0.6 && 0.019 < wb < 0.025},
   {{Om, Om /. s}, {wb, wb /. s}}];
profL[caso_, "wb", v_, s_] := First@Quiet@FindMinimum[
   {chi2Full[Om, h, -1, 0, v, caso], 0.1 < Om < 0.6 && 0.5 < h < 0.9},
   {{Om, Om /. s}, {h, h /. s}}];

(* w0waCDM profiles (only BAO+CMB): I move the other four, with w0+wa<0.        *)
profW["Om", v_, s_] := First@Quiet@FindMinimum[
   {chi2Full[v, h, w0, wa, wb, "cmb"],
    0.5 < h < 0.9 && -3 < w0 < 1 && -3 < wa < 2 && w0 + wa < 0 && 0.019 < wb < 0.025},
   {{h, h /. s}, {w0, w0 /. s}, {wa, wa /. s}, {wb, wb /. s}}];
profW["h", v_, s_] := First@Quiet@FindMinimum[
   {chi2Full[Om, v, w0, wa, wb, "cmb"],
    0.1 < Om < 0.6 && -3 < w0 < 1 && -3 < wa < 2 && w0 + wa < 0 && 0.019 < wb < 0.025},
   {{Om, Om /. s}, {w0, w0 /. s}, {wa, wa /. s}, {wb, wb /. s}}];
profW["w0", v_, s_] := First@Quiet@FindMinimum[
   {chi2Full[Om, h, v, wa, wb, "cmb"],
    0.1 < Om < 0.6 && 0.5 < h < 0.9 && -3 < wa < 2 && v + wa < 0 && 0.019 < wb < 0.025},
   {{Om, Om /. s}, {h, h /. s}, {wa, wa /. s}, {wb, wb /. s}}];
profW["wa", v_, s_] := First@Quiet@FindMinimum[
   {chi2Full[Om, h, w0, v, wb, "cmb"],
    0.1 < Om < 0.6 && 0.5 < h < 0.9 && -3 < w0 < 1 && w0 + v < 0 && 0.019 < wb < 0.025},
   {{Om, Om /. s}, {h, h /. s}, {w0, w0 /. s}, {wb, wb /. s}}];
profW["wb", v_, s_] := First@Quiet@FindMinimum[
   {chi2Full[Om, h, w0, wa, v, "cmb"],
    0.1 < Om < 0.6 && 0.5 < h < 0.9 && -3 < w0 < 1 && -3 < wa < 2 && w0 + wa < 0},
   {{Om, Om /. s}, {h, h /. s}, {w0, w0 /. s}, {wa, wa /. s}}];

(* BAO profiles (LCDM): Om uses the analytic A (chi2BAO); h*rd uses the explicit A. *)
profBAOom[v_]      := chi2BAO[v, -1, 0];
profBAOhrd[v_, s_] := First@Quiet@FindMinimum[
   {chi2BAOA[Om, cc/(100 v)], 0.1 < Om < 0.6}, {{Om, Om /. s}}];

(* Profiles for BAO+CMB+SNe (chi2CS): LCDM (Om,H0,wb) and w0waCDM (all five).
   They are copies of profL/profW replacing chi2Full[..,"cmb"] with chi2CS, to be
   able to measure the 1D sigma of each parameter also in the case with SNe.   *)
profLCS["Om", v_, s_] := First@Quiet@FindMinimum[
   {chi2CS[v, h, -1, 0, wb], 0.5 < h < 0.9 && 0.021 < wb < 0.024},
   {{h, h /. s}, {wb, wb /. s}}];
profLCS["h", v_, s_] := First@Quiet@FindMinimum[
   {chi2CS[Om, v, -1, 0, wb], 0.1 < Om < 0.6 && 0.021 < wb < 0.024},
   {{Om, Om /. s}, {wb, wb /. s}}];
profLCS["wb", v_, s_] := First@Quiet@FindMinimum[
   {chi2CS[Om, h, -1, 0, v], 0.1 < Om < 0.6 && 0.5 < h < 0.9},
   {{Om, Om /. s}, {h, h /. s}}];
profWCS["Om", v_, s_] := First@Quiet@FindMinimum[
   {chi2CS[v, h, w0, wa, wb],
    0.5 < h < 0.9 && -3 < w0 < 1 && -3 < wa < 2 && w0 + wa < 0 && 0.021 < wb < 0.024},
   {{h, h /. s}, {w0, w0 /. s}, {wa, wa /. s}, {wb, wb /. s}}];
profWCS["h", v_, s_] := First@Quiet@FindMinimum[
   {chi2CS[Om, v, w0, wa, wb],
    0.1 < Om < 0.6 && -3 < w0 < 1 && -3 < wa < 2 && w0 + wa < 0 && 0.021 < wb < 0.024},
   {{Om, Om /. s}, {w0, w0 /. s}, {wa, wa /. s}, {wb, wb /. s}}];
profWCS["w0", v_, s_] := First@Quiet@FindMinimum[
   {chi2CS[Om, h, v, wa, wb],
    0.1 < Om < 0.6 && 0.5 < h < 0.9 && -3 < wa < 2 && v + wa < 0 && 0.021 < wb < 0.024},
   {{Om, Om /. s}, {h, h /. s}, {wa, wa /. s}, {wb, wb /. s}}];
profWCS["wa", v_, s_] := First@Quiet@FindMinimum[
   {chi2CS[Om, h, w0, v, wb],
    0.1 < Om < 0.6 && 0.5 < h < 0.9 && -3 < w0 < 1 && w0 + v < 0 && 0.021 < wb < 0.024},
   {{Om, Om /. s}, {h, h /. s}, {w0, w0 /. s}, {wb, wb /. s}}];
profWCS["wb", v_, s_] := First@Quiet@FindMinimum[
   {chi2CS[Om, h, w0, wa, v],
    0.1 < Om < 0.6 && 0.5 < h < 0.9 && -3 < w0 < 1 && -3 < wa < 2 && w0 + wa < 0},
   {{Om, Om /. s}, {h, h /. s}, {w0, w0 /. s}, {wa, wa /. s}}];

(* From profile to sigma by curvature at 2 points: if at a distance 'del' from the
   optimum the profile rises by Dchi^2 = dp (right) and dm (left), the parabola
   Dchi^2=(x/sigma)^2 gives  sigma = del*Sqrt[2/(dp+dm)] (with del=sigma it would
   give dp=dm=1). I take 'del' of the order of the expected sigma; the ratio
   corrects the rest. Only 2 evaluations of the profile per parameter.
   'c0' = minimum chi^2 = First[fit].                                          *)
sigCurv[prof_, best_?NumericQ, c0_?NumericQ, del_?NumericQ] := Module[{dp, dm},(*best is the value of the parameter at which the
minimum of chi^2 is reached (for example h=...) and c_0 is the value of chi^2 at that minimum.*)
   dp = Max[0, prof[best + del] - c0];
   dm = Max[0, prof[best - del] - c0];
   del Sqrt[2/(dp + dm)]]; (*We could use the same algorithm as for H_0, but here there are many parameters for which sigma has to be computed
   and in the previous one we iterated until finding a point above and below each deltachi=1 on the left and right (that is too costly). Instead,
   we know that
   $$\mathcal{L}(x) \propto \exp\left( -\frac{\chi^2(x)}{2} \right)$$
   On the other hand, by the Central Limit Theorem,
   we know that if we have enough data, the probability of the parameter $x$ must follow a Gaussian (or Normal) distribution centred
   on $x_{\text{best}}$ and with a standard deviation $\sigma$:
   $$\mathcal{L}(x) \propto \exp\left( -\frac{(x - x_{\text{best}})^2}{2\sigma^2} \right)$$
   Since both expressions describe exactly the same probability, their exponents must be equal, hence:
   $$\frac{\Delta\chi^2(x)}{2} = \frac{(x - x_{\text{best}})^2}{2\sigma^2}$$. Doing a 2nd-order Taylor expansion around the minimum:
   $$\chi^2(x) \approx \chi^2(x_{\text{best}}) + \left. \frac{d\chi^2}{dx} \right\vert{}_{x_{\text{best}}} (x - x_{\text{best}}) +
   \frac{1}{2} \left. \frac{d^2\chi^2}{dx^2} \right\vert{}_{x_{\text{best}}} (x - x_{\text{best}})^2 $$. Moving the order-0 term to the left,
   since $$\Delta\chi^2(x) = \chi^2(x) - \chi^2(x_{\text{best}})$$ and the first derivative vanishes at the minimum, we get:
   $$\Delta\chi^2(x) \approx \frac{1}{2} \left. \frac{d^2\chi^2}{dx^2} \right\vert{}_{x_{\text{best}}} (x - x_{\text{best}})^2$$. So we get:
   $$\frac{1}{2} \frac{d^2\chi^2}{dx^2} = \frac{1}{\sigma^2}$$. Setting Deltachi^2=1 and solving, we obtain:
   $$\sigma = \sqrt{ \frac{2}{\frac{d^2\chi^2}{dx^2}} }$$. Now, in numerical analysis, the second derivative of a function $\chi^2(x)$ at the central point $x$
   can be approximated discretely as:
   $$\frac{d^2\chi^2}{dx^2} \approx \frac{\chi^2(x+\delta) - 2\chi^2(x) + \chi^2(x-\delta)}{\delta^2}$$. Now:
   $$\text{dp} + \text{dm} = \chi^2(x+\delta) - 2\chi^2(x) + \chi^2(x-\delta)$$. Therefore,
   $$\sigma = \sqrt{ \frac{2}{\frac{\text{dp} + \text{dm}}{\delta^2}} } = \delta \sqrt{\frac{2}{\text{dp} + \text{dm}}}$$. It is enough to take delta
   small enough.*)

(* --- Sigma of each parameter (best = fitted value, c0 = First[fit],
       del of the order of the expected sigma, only as the step for the curvature). --- *)
sBAOom  = sigCurv[profBAOom, Om /. fitBAOl[[2]], First[fitBAOl], 0.008];
sBAOhrd = sigCurv[profBAOhrd[#, fitBAOl[[2]]] &,
                  cc/(100 Aopt[Om /. fitBAOl[[2]], -1, 0]), First[fitBAOl], 0.7];
sBBNom  = sigCurv[profL["bbn", "Om", #, fitBBNl[[2]]] &, Om /. fitBBNl[[2]], First[fitBBNl], 0.008];
sBBNh   = sigCurv[profL["bbn", "h",  #, fitBBNl[[2]]] &, h  /. fitBBNl[[2]], First[fitBBNl], 0.006];
sBBNwb  = sigCurv[profL["bbn", "wb", #, fitBBNl[[2]]] &, wb /. fitBBNl[[2]], First[fitBBNl], 0.0005];
sCMBom  = sigCurv[profL["cmb", "Om", #, fitCMBl[[2]]] &, Om /. fitCMBl[[2]], First[fitCMBl], 0.004];
sCMBh   = sigCurv[profL["cmb", "h",  #, fitCMBl[[2]]] &, h  /. fitCMBl[[2]], First[fitCMBl], 0.004];
sCMBwb  = sigCurv[profL["cmb", "wb", #, fitCMBl[[2]]] &, wb /. fitCMBl[[2]], First[fitCMBl], 0.00012];
sWom = sigCurv[profW["Om", #, fitCMBw[[2]]] &, Om /. fitCMBw[[2]], First[fitCMBw], 0.024];
sWh  = sigCurv[profW["h",  #, fitCMBw[[2]]] &, h  /. fitCMBw[[2]], First[fitCMBw], 0.021];
sWw0 = sigCurv[profW["w0", #, fitCMBw[[2]]] &, w0 /. fitCMBw[[2]], First[fitCMBw], 0.24];
sWwa = sigCurv[profW["wa", #, fitCMBw[[2]]] &, wa /. fitCMBw[[2]], First[fitCMBw], 0.7];
sWwb = sigCurv[profW["wb", #, fitCMBw[[2]]] &, wb /. fitCMBw[[2]], First[fitCMBw], 0.00014];

(* One-dimensional sigmas of the BAO+CMB+SNe case: LCDM (Om,H0,wb) and w0waCDM (all five). *)
sCSom  = sigCurv[profLCS["Om", #, fitCSl[[2]]] &, Om /. fitCSl[[2]], First[fitCSl], 0.004];
sCSh   = sigCurv[profLCS["h",  #, fitCSl[[2]]] &, h  /. fitCSl[[2]], First[fitCSl], 0.004];
sCSwb  = sigCurv[profLCS["wb", #, fitCSl[[2]]] &, wb /. fitCSl[[2]], First[fitCSl], 0.00012];
sWCSom = sigCurv[profWCS["Om", #, fitCSw[[2]]] &, Om /. fitCSw[[2]], First[fitCSw], 0.006];
sWCSh  = sigCurv[profWCS["h",  #, fitCSw[[2]]] &, h  /. fitCSw[[2]], First[fitCSw], 0.006];
sWCSw0 = sigCurv[profWCS["w0", #, fitCSw[[2]]] &, w0 /. fitCSw[[2]], First[fitCSw], 0.06];
sWCSwa = sigCurv[profWCS["wa", #, fitCSw[[2]]] &, wa /. fitCSw[[2]], First[fitCSw], 0.2];
sWCSwb = sigCurv[profWCS["wb", #, fitCSw[[2]]] &, wb /. fitCSw[[2]], First[fitCSw], 0.00014];

(* Formats a cell  "value +- sigma"  with the same number format in both. *)
fmtPM[v_, e_, f_] := ToString@NumberForm[v, f] <> " \[PlusMinus] " <> ToString@NumberForm[e, f];

(* The "DESI DR2" column also carries the +-sigma published by DESI, and the
   decimals of "mi ajuste" have been matched to those of DESI so that one can
   compare at a glance. Where the DESI sigmas come from:
     - BAO (LCDM) and BAO+BBN (LCDM): eq.17 and eq.19 of the text.
     - BAO+CMB (LCDM): Table 5 (the one I used to call "Table 3"; in the
       cells it is already corrected to Table 5).
     - BAO+CMB (w0waCDM): eq.24 (w0, wa, Om) and its H0.
     - wb: in BBN it is the BBN prior (0.02218+-0.00055); in CMB the DESI column
       carries the Planck value (0.02236+-0.00015). Note that my compressed vector
       now uses the wCDM row (0.02239); I leave the DESI reference at 0.02236.
   Caveat 1: the H0 of DESI+CMB in w0waCDM is asymmetric (about +1.6/-2.1); I put
   it as approximately +-1.9.
   Caveat 2: DESI compares with its full CMB (Planck likelihood), not with
   compressed distance priors like mine; that is why I compare my BAO+CMB
   (compressed) against the full DESI+CMB, which is the right thing to do (see
   Table B).
   Confirm the exact sigmas in Table 5 before submitting; the central values
   already matched.                                                           *)
(* --- chi^2_min/dof (reduced chi^2) per fit, for Tables A1 and A2 --
   I add to A1 (LCDM) and A2 (w0waCDM) a column with the reduced chi^2 =
   chi^2_min/dof of each fit. It is a single cell per case, which spans its
   parameter rows with SpanFromAbove (like the "Caso / modelo" column), so it
   takes one line per cell. The chi^2_min is already computed in Block 3 (it is
   First[fit...]); only the dof is missing here.

     dof = (number of measurements entering the chi^2) - (number of free parameters).

   (1) Measurements per case (terms adding to chi^2_min; see Blocks 2 and 3):
        BAO          : 13         (= Length[dvec], the 13 DESI observables)
        BAO+BBN      : 13 + 1     (+1 = Gaussian prior on wb; it adds one term)
        BAO+CMB      : 13 + 3     (+3 = distance priors R, l_A, wb)
        BAO+CMB+SNe  : 13 + 3 + 7 (+7 = the 7 Pantheon+ 1/E(z) nodes)
      Convention: I count as a "measurement" each Gaussian term of the chi^2,
      including the BBN/CMB priors, which is consistent with chi^2_min already
      including them. If one prefers not to count the priors as data, it is enough
      to change the nDat* summands below (the only convention here).
   (2) Free parameters (the ones NMinimize varies in Block 3):
        LCDM   : BAO -> {Om, A} = 2  (A is marginalised in chi2BAO, but it uses up one dof)
                 BAO+BBN / BAO+CMB / BAO+CMB+SNe -> {Om, h, wb} = 3
        w0waCDM: (in A2 only BAO+CMB and BAO+CMB+SNe) -> {Om, h, w0, wa, wb} = 5 . *)
nBAO    = Length[dvec];              (* = 13 BAO observables.                     *)
nDatBAO = nBAO;                      (* BAO         : 13.                         *)
nDatBBN = nBAO + 1;                  (* BAO+BBN     : 13 + 1 (prior wb).          *)
nDatCMB = nBAO + 3;                  (* BAO+CMB     : 13 + 3 (R, l_A, wb).        *)
nDatCS  = nBAO + 3 + 7;              (* BAO+CMB+SNe : 13 + 3 + 7 (7 SNe nodes).   *)
dofBAOl = nDatBAO - 2;               (* BAO      LCDM    : 13-2 = 11.             *)
dofBBNl = nDatBBN - 3;               (* BAO+BBN  LCDM    : 14-3 = 11.             *)
dofCMBl = nDatCMB - 3;               (* BAO+CMB  LCDM    : 16-3 = 13.             *)
dofCSl  = nDatCS  - 3;               (* +SNe     LCDM    : 23-3 = 20.             *)
dofCMBw = nDatCMB - 5;               (* BAO+CMB  w0waCDM : 16-5 = 11.             *)
dofCSw  = nDatCS  - 5;               (* +SNe     w0waCDM : 23-5 = 18.             *)

(* Cell of the new column: "chi^2_min / dof = reduced chi^2", in one line. *)
redCell[chi2_, dof_] := Row[{NumberForm[chi2, {4,1}], "/", dof, " = ",
     NumberForm[chi2/dof, {4,2}]}];

tablaLCDM = {
  {"BAO (\[CapitalLambda]CDM)",     redCell[First[fitBAOl], dofBAOl],   (* Chi^2/dof of the fit. *)
      "\[CapitalOmega]m",
      fmtPM[Om /. fitBAOl[[2]], sBAOom, {5,4}],                         "0.2975 \[PlusMinus] 0.0086"},
  {SpanFromAbove,                    SpanFromAbove,                     (* The chi^2/dof cell spans the case. *)
      "h\[CenterDot]rd",
      fmtPM[cc/(100 Aopt[Om /. fitBAOl[[2]], -1, 0]), sBAOhrd, {5,2}],  "101.54 \[PlusMinus] 0.73"},
  {"BAO+BBN (\[CapitalLambda]CDM)", redCell[First[fitBBNl], dofBBNl],   
      "\[CapitalOmega]m",
      fmtPM[Om /. fitBBNl[[2]], sBBNom, {5,4}],                         "0.2977 \[PlusMinus] 0.0086"},
  {SpanFromAbove,                    SpanFromAbove,                     
      "H0",
      fmtPM[100 (h /. fitBBNl[[2]]), 100 sBBNh, {4,2}],                 "68.51 \[PlusMinus] 0.58"},
  {SpanFromAbove,                    SpanFromAbove,                     
      "\[Omega]b",
      fmtPM[wb /. fitBBNl[[2]], sBBNwb, {6,5}],                         "0.02218 \[PlusMinus] 0.00055"},
  {"BAO+CMB (\[CapitalLambda]CDM)", redCell[First[fitCMBl], dofCMBl],   
      "\[CapitalOmega]m",
      fmtPM[Om /. fitCMBl[[2]], sCMBom, {4,3}],                         "0.301 \[PlusMinus] 0.004"},
  {SpanFromAbove,                    SpanFromAbove,                     
      "H0",
      fmtPM[100 (h /. fitCMBl[[2]]), 100 sCMBh, {4,2}],                 "68.44 \[PlusMinus] 0.30"},
  {SpanFromAbove,                    SpanFromAbove,                     
      "\[Omega]b",
      fmtPM[wb /. fitCMBl[[2]], sCMBwb, {6,5}],                         "0.02236 \[PlusMinus] 0.00015"},
  (* BAO+CMB+SNe in LCDM (the most complete case). DESI does not tabulate LCDM with
     BAO+CMB+Pantheon+, so the DESI cells for Om and H0 only indicate that; in
     LCDM the SNe barely move Om and H0 with respect to BAO+CMB. *)
  {"BAO+CMB+SNe (\[CapitalLambda]CDM)", redCell[First[fitCSl], dofCSl], 
      "\[CapitalOmega]m",
      fmtPM[Om /. fitCSl[[2]], sCSom, {4,3}],                           "DESI no lo tabula"},
  {SpanFromAbove,                    SpanFromAbove,                     
      "H0",
      fmtPM[100 (h /. fitCSl[[2]]), 100 sCSh, {4,2}],                   "DESI no lo tabula"},
  {SpanFromAbove,                    SpanFromAbove,                     
      "\[Omega]b",
      fmtPM[wb /. fitCSl[[2]], sCSwb, {6,5}],                           "0.02236 \[PlusMinus] 0.00015"}};
Print["\nTABLA A1 : parametros \[CapitalLambda]CDM  (mi ajuste  vs  DESI DR2)"];
estiloTabla[tablaLCDM, {"Caso / modelo",
   Row[{Subsuperscript["\[Chi]", "min", "2"], "/dof"}],   
   "parametro", "mi ajuste (valor \[PlusMinus] \[Sigma])", "DESI DR2"}]

tablaW0WA = {
  {"BAO+CMB (w0waCDM)", redCell[First[fitCMBw], dofCMBw],   (* chi^2/dof of the fit. *)
      "\[CapitalOmega]m",
      fmtPM[Om /. fitCMBw[[2]], sWom, {4,3}],                           "0.353 \[PlusMinus] 0.021"},
  {SpanFromAbove,        SpanFromAbove,                     (* The chi^2/dof cell spans the case. *)
      "H0",
      fmtPM[100 (h /. fitCMBw[[2]]), 100 sWh, {3,1}],                   "63.6 \[PlusMinus] 1.9"},
  {SpanFromAbove,        SpanFromAbove,                     
      "w0",
      fmtPM[w0 /. fitCMBw[[2]], sWw0, {3,2}],                           "-0.42 \[PlusMinus] 0.21"},
  {SpanFromAbove,        SpanFromAbove,                     
      "wa",
      fmtPM[wa /. fitCMBw[[2]], sWwa, {3,2}],                           "-1.75 \[PlusMinus] 0.58"},
  {SpanFromAbove,        SpanFromAbove,                     
      "\[Omega]b",
      fmtPM[wb /. fitCMBw[[2]], sWwb, {6,5}],                           "0.02236 \[PlusMinus] 0.00015"},
  (* BAO+CMB+SNe in w0waCDM, with the DESI+CMB+Pantheon+ values
     (DESI DR2, 2503.14738; PP = Pantheon+). The decimals of "mi ajuste" are
     matched to those of DESI. The wa of DESI is asymmetric (+0.22/-0.19). *)
  {"BAO+CMB+SNe (w0waCDM)", redCell[First[fitCSw], dofCSw],   (* chi^2/dof of the fit. *)
      "\[CapitalOmega]m",
      fmtPM[Om /. fitCSw[[2]], sWCSom, {5,4}],                          "0.3114 \[PlusMinus] 0.0057"},
  {SpanFromAbove,        SpanFromAbove,                     
      "H0",
      fmtPM[100 (h /. fitCSw[[2]]), 100 sWCSh, {4,2}],                  "67.51 \[PlusMinus] 0.59"},
  {SpanFromAbove,        SpanFromAbove,                     
      "w0",
      fmtPM[w0 /. fitCSw[[2]], sWCSw0, {4,3}],                          "-0.838 \[PlusMinus] 0.055"},
  {SpanFromAbove,        SpanFromAbove,                     
      "wa",
      fmtPM[wa /. fitCSw[[2]], sWCSwa, {3,2}],                          "-0.62 +0.22/-0.19"},
  {SpanFromAbove,        SpanFromAbove,                     
      "\[Omega]b",
      fmtPM[wb /. fitCSw[[2]], sWCSwb, {6,5}],                          "0.02236 \[PlusMinus] 0.00015"}};
Print["\nTABLA A2 : parametros w0waCDM = SOLO BAO+CMB  (mi ajuste  vs  DESI DR2)"];
estiloTabla[tablaW0WA, {"Caso / modelo",
   Row[{Subsuperscript["\[Chi]", "min", "2"], "/dof"}],   
   "parametro", "mi ajuste (valor \[PlusMinus] \[Sigma])", "DESI DR2"}]

(* Table B: significance of w0waCDM against LCDM ---------------------- *)
tablaSig = {
  {"BAO",     ToString@NumberForm[DchiBAO, {3,1}],
              ToString@NumberForm[sigDe[DchiBAO], {2,1}] <> "\[Sigma]",
              "~1.7\[Sigma]"},
  {"BAO+BBN", ToString@NumberForm[DchiBBN, {3,1}],
              ToString@NumberForm[sigDe[DchiBBN], {2,1}] <> "\[Sigma]",
              "DESI no lo tabula suelto"},
  {"BAO+CMB", ToString@NumberForm[DchiCMB, {3,1}],
              ToString@NumberForm[sigDe[DchiCMB], {2,1}] <> "\[Sigma]",
              "3.1\[Sigma]"},
  {"BAO+CMB+SNe", ToString@NumberForm[DchiCS, {3,1}],
              ToString@NumberForm[sigDe[DchiCS], {2,1}] <> "\[Sigma]",
              "2.8\[Sigma]"}};
Print["\nTABLA B : significancia frente a \[CapitalLambda]CDM"];
estiloTabla[tablaSig,
   {"Caso", Row[{"\[CapitalDelta]", Superscript["\[Chi]", "2"], " MAP"}],
    "significancia", "DESI DR2"}]


(* =====================================================================
   Block 6 : how to compare with DESI (quick check)
   ---------------------------------------------------------------------
   1) LCDM: Om ~ 0.297 and (with BBN) H0 ~ 68.5, see eqs. 17 and 19.
   2) w0waCDM+CMB: the best (w0,wa) is in the quadrant w0>-1, wa<0, with Om~0.35,
      see eq. 24 and Table V.
   3) Significance: Dchi^2_MAP and sigmas, comparable to Table VI. My compressed
      CMB gives about 2.3 sigma; the 3.1 sigma of DESI uses the full CMB (it adds
      lensing and information that the three-number compression does not capture).
   4) LCDM in the (w0,wa) plane: inside the 68-95% with BAO (compatible),
      but outside the 95% with BAO+CMB (preference for dynamical dark
      energy), which is why there I also draw the 3 sigma contour.
   5) w0waCDM+CMB+SNe: adding Pantheon+ (compressed) the fit sticks closely to
      DESI+CMB+Pantheon+ (w0~-0.84, wa~-0.6, Om~0.31, H0~67.5); DESI gives 2.8
      sigma and my compressed version about 2.2 sigma (see Tables A2 and B).
   --------------------------------------------------------------------- *)


(* =====================================================================
   Block 7 : the supernovae (how they are included: compressed Pantheon+)
   ---------------------------------------------------------------------
   They are already included (the BAO+CMB+SNe case). Instead of the full vector of
   distance moduli mu_obs(z_i) of a catalogue and its covariance, I use the
   compression of Pantheon+ (Table II of 2408.17318): 1/E(z) is fitted as a
   spline to the full Pantheon+ and the script-M nuisance (mix of M_B and H0)
   is marginalised, giving 7 measurements of 1/E(z) with their 7x7 covariance.
   The SNe chi^2 becomes
        chi^2_SNe(Om,w0,wa) = (g - eSN)^T CinvSN (g - eSN),  g_i = 1/E(z_i),
   which only depends on (Om,w0,wa) (the script-M is already marginalised, like
   the A of BAO) and is added to BAO+CMB, giving chi2CS (Block 2). The input
   data (zSN, eSN, sSN, rhoSN) are in Block 0.
   Consistency: my fit (w0~-0.86, wa~-0.50, Om~0.31, H0~67.6, about 2.2 sigma)
   nearly matches the published DESI+CMB+Pantheon+. For Union3 or DES-SN5YR
   their own analogous compressed table would suffice.
   ===================================================================== *)


(* ==== Block 8 : export of the plots to PDF ==== *)

(* 1. We tell Mathematica to work in the same folder where this file is saved. *)
SetDirectory[NotebookDirectory[]];

(* 2. We export the 7 figures to high-quality vector PDF format. *)
Export["A_Contornos_BAO.pdf", figBAO];
Export["B_Contornos_BAO_BBN.pdf", figBBNcont];
Export["C_Contornos_BAO_CMB.pdf", figCMB];
Export["D_Contornos_BAO_CMB_SNe.pdf", figCS];
Export["E_Tension_Hubble_H0.pdf", figBBN];
Export["F_Restriccion_OmegaM_Verosimilitud.pdf", figOm];
Export["G_Contornos_Solape_w0wa.pdf", figSolape];

(* 3. Confirmation message, to know that everything went well. *)
Print["\[DownExclamation]Exito! Las 7 graficas del ajuste DESI se han guardado en la carpeta: ", Directory[]];
