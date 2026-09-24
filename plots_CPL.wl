(* ::Package:: *)

(* =====================================================================
   plots_CPL.wl   (Wolfram Language / Mathematica)
   ---------------------------------------------------------------------
   Plots of the expansion history comparing three cosmologies:
        - fiducial LCDM            (Om=0.31, w=-1).
        - w0waCDM  BAO+CMB         (best fit).
        - w0waCDM  BAO+CMB+SNe     (best fit).
   with the CPL parametrisation  w(a)=w0+wa(1-a).

   File separate from the fitting one (DESI_fit.wl): that one was already
   loaded with calculations; here we only plot with the already fitted parameters.

   Figures:
     A) E(z)   (3 lines)      C) a(t)   (3 lines)                               E) fractions Omega_i.
     B) q(z)   (3 lines)      D) densities rho_i (CPL BAO+CMB and BAO+CMB+SNe)
   -> A,B,C in a single plot each; D and E are 2 plots each (one
      per w0waCDM model), because the components of two models do not fit
      legibly in the same panel.
   (The distances plot is kept only in LCDM, in your cosmological distances
    section; not included here since it does not separate the models: see Notes.)
   ===================================================================== *)


(* ==== BLOCK 1 : constants ==== *)
cc = 299792.458;          (* c [km/s]. *)
wr = 4.178*^-5;           (* omega_r = Omega_r h^2 (photons+neutrinos), as in the fit. *)


(* ==== BLOCK 2 : opimal parameters (change them here if your fit varies) ====
   {Om, w0, wa, h}. Omega_r of each model = wr/h^2 (the radiation today as a
   fraction depends on h). LCDM goes with w0=-1, wa=0 (constant dark energy). *)
parLCDM = {0.31,   -1.0,     0.0,     0.70};   (* Fiducial from the plotting code. *)
parCMB  = {0.3488, -0.4723, -1.5618, 0.6399};  (* w0waCDM  BAO+CMB.       *)
parCS   = {0.3113, -0.8554, -0.5001, 0.6758};  (* w0waCDM  BAO+CMB+SNe.   *)
orrOf[h_] := wr/h^2;                            (* Omega_r today = wr/h^2. *)


(* ==== BLOCK 3 : the CPL model (E, w, q and DE density) ====
   CPL:  w(a)=w0+wa(1-a)  ->  in z:  w(z)=w0+wa z/(1+z).
   Dark energy density (normalised to today), integrating the continuity
   eq. of a fluid with w(a):
        rho_DE(z)/rho_DE0 = (1+z)^{3(1+w0+wa)} Exp[-3 wa z/(1+z)].
   (for LCDM, w0=-1,wa=0 -> it equals 1: constant Lambda).
   E(z)^2 = Or(1+z)^4 + Om(1+z)^3 + (1-Om-Or) rho_DE(z)/rho_DE0.               *)
fDE[z_, w0_, wa_] := (1 + z)^(3 (1 + w0 + wa)) Exp[-3 wa z/(1 + z)]; (* rho_DE/rho_DE0. *)
EzC[z_, Om_, w0_, wa_, orr_] :=
   Sqrt[orr (1 + z)^4 + Om (1 + z)^3 + (1 - Om - orr) fDE[z, w0, wa]];
wC[z_, w0_, wa_] := w0 + wa z/(1 + z);                              (* w(z). *)
(* q(z) = sum_i Omega_i(z) (1+3 w_i)/2 ; radiation w=1/3, matter w=0, DE w(z). *)
qC[z_, Om_, w0_, wa_, orr_] := Module[{e2, fr, fm, fd},
   e2 = orr (1 + z)^4 + Om (1 + z)^3 + (1 - Om - orr) fDE[z, w0, wa];
   fr = orr (1 + z)^4/e2; fm = Om (1 + z)^3/e2; fd = (1 - Om - orr) fDE[z, w0, wa]/e2;
   fr + fm/2 + fd (1 + 3 wC[z, w0, wa])/2];


(* ==== BLOCK 3b : format+manual legend (must be evaluated) ====
   The legend swatch is a real line (not a glyph),
   so it never comes out broken, and it all goes in one Row so it does not split. *)
colLCDM = RGBColor[0.106, 0.286, 0.396];   (* navy blue. *)
colCMB  = RGBColor[0.820, 0.150, 0.150];   (* red.         *)
colCS   = RGBColor[0.150, 0.550, 0.300];   (* green.       *)
colMod  = {colLCDM, colCMB, colCS};
legMod  = {"\[CapitalLambda]CDM", "BAO+CMB", "BAO+CMB+SNe"};

(* colours of the components (densities/fractions). *)
colRad = RGBColor[0.894, 0.631, 0.106];    (* radiation. *)
colMat = RGBColor[0.106, 0.286, 0.396];    (* matter.    *)
colDE  = RGBColor[0.706, 0.271, 0.122];    (* dark energy. *)

sH0  = Subscript["H", "0"];
sOmi = Subscript["\[CapitalOmega]", "i"];
sRhoi = Subscript["\[Rho]", "i"];
sRhoc = Subscript["\[Rho]", "c,0"];

muestra[col_, dash_:Nothing] := Graphics[
   {col, AbsoluteThickness[2.4], dash, Line[{{0, 0}, {1, 0}}]},
   ImageSize -> {30, 10}, AspectRatio -> Full, PlotRange -> {{0, 1}, {-1, 1}},
   PlotRangePadding -> None, ImagePadding -> None, BaselinePosition -> Center,
   Frame -> False, Axes -> False, GridLines -> None, Background -> None];
itemLey[col_, lab_, dash_ : Nothing] := Row[{muestra[col, dash], Spacer[5], Style[lab, Black, 12]}];
filaLey[cols_, labs_] := Row[Riffle[MapThread[itemLey, {cols, labs}], Spacer[22]]];
(* [NEW] variant with a list of DASH patterns (one per curve) to tell them apart in B&W *)
filaLey[cols_, labs_, dashes_] := Row[Riffle[MapThread[itemLey, {cols, labs, dashes}], Spacer[22]]];

(* -- Ticks on the logarithmic axes, all as 10^k ----------------------
   By default Mathematica mixes on the same logarithmic axis the decimal
   notation (0.001, 0.1, 10, 1000) with the scientific one (10^-5, 10^5, ...).
   Here the ticks are built by hand so that all of them come out as 10^k.

     marcasLog[kmin, kmax]        -> label 10^k at every decade
     marcasLog[kmin, kmax, paso]  -> label only every 'paso' decades; the
                                     ones in between stay as short ticks
                                     WITHOUT label (so text does not pile up)
     marcasLogMudas[...]          -> the same ticks without labels, for the
                                     top and right sides of the frame.        *)
marcasLog[kmin_Integer, kmax_Integer, paso_Integer : 1] :=
  Module[{conEtiq, sinEtiq},
   conEtiq = Range[kmin, kmax, paso];                  (* labelled decades. *)
   sinEtiq = Complement[Range[kmin, kmax], conEtiq];   (* intermediate decades. *)
   Join[
     Table[{10.^k, Superscript[10, k], {0.010, 0}}, {k, conEtiq}],
     Table[{10.^k, "", {0.005, 0}}, {k, sinEtiq}]]];

marcasLogMudas[kmin_Integer, kmax_Integer, paso_Integer : 1] :=
  ({#[[1]], "", #[[3]]} &) /@ marcasLog[kmin, kmax, paso];

(* -- ticks of a linear axis (the Omega_i one, from 0 to 1). They are also
      written by hand only to be able to repeat them without label on the
      right side, as the default format does. Labels 0.0, 0.2, ..., 1.0.       *)
marcasLin[min_, max_, mayor_, menor_] :=
  Module[{conEtiq, sinEtiq},
   conEtiq = Range[min, max, mayor];
   sinEtiq = Complement[Range[min, max, menor], conEtiq];
   Join[
     Table[{x, NumberForm[N[x], {3, 1}], {0.010, 0}}, {x, conEtiq}],
     Table[{x, "", {0.005, 0}}, {x, sinEtiq}]]];

marcasLinMudas[min_, max_, mayor_, menor_] :=
  ({#[[1]], "", #[[3]]} &) /@ marcasLin[min, max, mayor, menor];

SetOptions[{Plot, LogLogPlot, LogLinearPlot, ParametricPlot},
   Frame -> True, Axes -> False, GridLines -> Automatic,
   GridLinesStyle -> Directive[GrayLevel[0.85], AbsoluteThickness[0.4]],
   FrameStyle -> Directive[Black, AbsoluteThickness[1]],
   LabelStyle -> Directive[Black, 12], Background -> White,
   ImageSize -> 500, AspectRatio -> 0.72, PlotRangePadding -> Scaled[0.015],
   ImagePadding -> {{Automatic, 22}, {Automatic, Automatic}}];

(* Dash pattern per model: solid/dashed/dotted (B&W/colour-blind).
   dash3 is also reused by the components of figDens/figFrac. *)
dash3 = {Dashing[None], AbsoluteDashing[{9, 5}], AbsoluteDashing[{1, 4}]};
(* est3 combines colour + thickness + dash pattern. *)
est3 = MapThread[Directive[#1, AbsoluteThickness[2], #2] &, {colMod, dash3}];
(* dash-dot style (dotted-dashed) for the equality vertical lines *)
verticalEq = Directive[GrayLevel[0.35], AbsoluteDashing[{6, 3, 1, 3}], AbsoluteThickness[1.3]];


(* ==== BLOCK 4 -> FIGURE A : E(z) (3 models) ==== *)
figEz = Legended[
  Plot[Evaluate@{
     EzC[z, parLCDM[[1]], parLCDM[[2]], parLCDM[[3]], orrOf[parLCDM[[4]]]],
     EzC[z, parCMB[[1]],  parCMB[[2]],  parCMB[[3]],  orrOf[parCMB[[4]]]],
     EzC[z, parCS[[1]],   parCS[[2]],   parCS[[3]],   orrOf[parCS[[4]]]]},
    {z, 0, 3}, PlotStyle -> est3,
    FrameLabel -> {"z", Row[{"E(z) = H/", sH0}]}],
  Placed[filaLey[colMod, legMod, dash3], Below]]


(* ==== BLOCK 5 -> FIGURE B : q(z) (3 models) ====
   q=0 marks the onset of acceleration; I draw that horizontal line. *)
zacc[p_] := z /. FindRoot[qC[z, p[[1]], p[[2]], p[[3]], orrOf[p[[4]]]] == 0, {z, 0.6}];
figQz = Legended[
  Plot[Evaluate@{
     qC[z, parLCDM[[1]], parLCDM[[2]], parLCDM[[3]], orrOf[parLCDM[[4]]]],
     qC[z, parCMB[[1]],  parCMB[[2]],  parCMB[[3]],  orrOf[parCMB[[4]]]],
     qC[z, parCS[[1]],   parCS[[2]],   parCS[[3]],   orrOf[parCS[[4]]]]},
    {z, 0, 3}, PlotStyle -> est3,
    FrameLabel -> {"z", "q(z)"},
    PlotRange -> {Automatic, {-0.6, Automatic}}, 
    Epilog -> {Directive[GrayLevel[0.4], Dashed, AbsoluteThickness[1.3]],
               Line[{{0, 0}, {3, 0}}]}],
  Placed[filaLey[colMod, legMod, dash3], Below]]


(* ==== BLOCK 6 -> FIGURE C : a(t) (3 models) ====
   t-t0 = (9.778/h) * Integral_1^a da'/(a' E(a')).  Each model with its own h
   (its H0), which is why the curves also diverge in time scale.              *)
tofa[a_?NumericQ, p_] := (9.778/p[[4]]) NIntegrate[
   1/(x EzC[1/x - 1, p[[1]], p[[2]], p[[3]], orrOf[p[[4]]]]), {x, 1, a}];
figAt = Legended[
  ParametricPlot[Evaluate@{
     {tofa[a, parLCDM], a}, {tofa[a, parCMB], a}, {tofa[a, parCS], a}},
    {a, 0.001, 3.7}, PlotStyle -> est3, PlotPoints -> 120, MaxRecursion -> 4,   (* a from 0.001 to 3.7 and fine sampling: draws the whole almost 
    vertical stretch of the Big Bang down to a=0 (at the bottom) and reaches t=20 (at the top), without being cut off. *)
    FrameLabel -> {Row[{"t - ", Subscript["t", "0"], " [Gyr]  (0 = hoy)"}], "a(t)"},
    PlotRange -> {{-15, 20}, {0, 3.6}},   (* Y axis up to 3.6: LCDM reaches a~3.5 at t=20. *)
    Epilog -> {{Directive[GrayLevel[0.4], Dotted, AbsoluteThickness[1.5]],
                Line[{{-15, 1}, {20, 1}}]},
               {Directive[GrayLevel[0.4], Dotted, AbsoluteThickness[1.5]],
                Line[{{0, 0}, {0, 3.6}}]}}],   (* Vertical line t=0 up to the new y maximum. *)
  Placed[filaLey[colMod, legMod, dash3], Below]]


(* ==== BLOCK 7 -> FIGURES D : densities rho_i(a)/rho_c,0 (CPL, per model) ====
   rho_rad = Or a^-4 ; rho_mat = Om a^-3 ; rho_DE = (1-Om-Or) fDE(z(a)) with
   a=1/(1+z). Now dark energy is not flat (CPL): it rises, peaks near
   today and decays (CPL extrapolation outside the data range).
   >>> As in the fractions, I mark with two dash-dot vertical lines the
       a values of density equality, with their value of a next to the axis. *)

figDens[p_, nombre_, xMinExp_, yMinExp_] := Module[
  {Om = p[[1]], w0 = p[[2]], wa = p[[3]], orr = orrOf[p[[4]]], aRM, aMDE, yTickStart},
  aRM  = orr/Om;                                                          (* radiation-matter. *)
  aMDE = aa /. FindRoot[Om aa^-3 == (1 - Om - orr) fDE[1/aa - 1, w0, wa], {aa, 0.7}]; (* matter-DE. *)
  
  (* We compute the start of the grid so that it always fits in multiples of 3. *)
  yTickStart = Floor[yMinExp / 3] * 3;
  
  Legended[
   LogLogPlot[{orr a^-4, Om a^-3, (1 - Om - orr) fDE[1/a - 1, w0, wa]}, {a, 10.^xMinExp, 3},
     PlotStyle -> {Directive[colRad, AbsoluteThickness[1.8], dash3[[1]]],
                   Directive[colMat, AbsoluteThickness[1.8], dash3[[2]]],
                   Directive[colDE,  AbsoluteThickness[1.8], dash3[[3]]]}, 
     PlotRange -> {{10.^xMinExp, 3}, {10.^yMinExp, 1*^20}},
     FrameLabel -> {"a", Row[{sRhoi, " / ", sRhoc}]},
     
     (* Uses yTickStart to keep the grid consistent *)
     FrameTicks -> {{marcasLog[yTickStart, 20, 3], marcasLogMudas[yTickStart, 20, 3]},
                    {marcasLog[xMinExp,  0],    marcasLogMudas[xMinExp,  0]}},
     GridLines -> {Join[Table[10.^k, {k, xMinExp, 0}],
                        {{aRM, verticalEq}, {aMDE, verticalEq}}],
                   Table[10.^k, {k, yTickStart, 18, 3}]},
                   
     Epilog -> {
        (* Intersection labels automatically adjusted to xMinExp *)
        Text[Framed[Style[Row[{"a\[TildeTilde]", ScientificForm[aRM, 2]}], 10, GrayLevel[0.25]],
              Background -> White, FrameStyle -> None, FrameMargins -> 1],
           Scaled[{(Log[10, aRM] - xMinExp)/(Log[10, 3] - xMinExp), 0.05}], {-1.6, 0}],
           
        Text[Framed[Style[Row[{"a\[TildeTilde]", NumberForm[aMDE, {3, 2}]}], 10, GrayLevel[0.25]],
              Background -> White, FrameStyle -> None, FrameMargins -> 1],
           Scaled[{(Log[10, aMDE] - xMinExp)/(Log[10, 3] - xMinExp), 0.06}], {1.8, 0}]}],
           
   Placed[filaLey[{colRad, colMat, colDE},
      {Row[{"Radiaci\[OAcute]n \[Proportional] ", Superscript["a", "-4"]}],
       Row[{"Materia \[Proportional] ", Superscript["a", "-3"]}],
       "DE (CPL)"}, dash3], Below]]];

(* Application of the limits for each model:
   - Plot 1 (CMB): X axis at -5, Y axis goes down to -14 (full radiation curve).
   - Plot 2 (CS):  X axis at -5, Y axis goes down to -5. *)
figDensCMB = figDens[parCMB, "BAO+CMB", -5, -14]
figDensCS  = figDens[parCS,  "BAO+CMB+SNe", -5, -5]


(* ==== BLOCK 8 -> fugures E : fractions Omega_i(z) (CPL, per model) ====
   Omega_i(z) = rho_i(z)/rho_total(z) ; they add up to 1. x axis = 1+z in log.
   >>> I mark with two dash-dot vertical lines the z of equality:
       - matter-radiation:   rho_m=rho_r  ->  1+z = Om/Or   (z ~ 3400)
       - matter-DE:          rho_m=rho_DE ->  solved with FindRoot (CPL is not
         a pure power law, so there is no closed formula; z ~ 0.3).
       Each vertical line carries next to the axis the VALUE OF z of the equality.
   The vertical lines go as GridLines (data coordinates, robust on a log axis)
   and the labels in Scaled (fraction = Log10[1+z]/7, with axis 1..1e7).       *)
figFrac[p_, nombre_] := Module[
  {Om = p[[1]], w0 = p[[2]], wa = p[[3]], orr = orrOf[p[[4]]], zMR, zMDE, xMR, xMDE},
  zMR  = Om/orr - 1;                                                   (* matter-radiation. *)
  zMDE = zz /. FindRoot[Om (1 + zz)^3 == (1 - Om - orr) fDE[zz, w0, wa], {zz, 0.3}]; (* matter-DE. *)
  xMR = 1 + zMR; xMDE = 1 + zMDE;
  Legended[
   LogLinearPlot[
     Evaluate@With[{e2 = orr u^4 + Om u^3 + (1 - Om - orr) fDE[u - 1, w0, wa]},
        {orr u^4/e2, Om u^3/e2, (1 - Om - orr) fDE[u - 1, w0, wa]/e2}],
     {u, 1, 1*^7},
     PlotStyle -> {Directive[colRad, AbsoluteThickness[1.8], dash3[[1]]],
                   Directive[colMat, AbsoluteThickness[1.8], dash3[[2]]],
                   Directive[colDE,  AbsoluteThickness[1.8], dash3[[3]]]}, (* Dash pattern per component (B&W). *)
     PlotRange -> {{1, 1*^7}, {-0.02, 1.02}},
     FrameLabel -> {"1+z", Row[{sOmi, "(z)"}]},
     (* X axis (log) ALL as 10^k: 10^0, 10^1, ..., 10^7, instead of
        1, 10, 100, 1000, 10^4, ... (which mixes decimal and scientific).
        The y axis is linear and stays as always: 0.0, 0.2, ..., 1.0.          *)
     FrameTicks -> {{marcasLin[0, 1, 1/5, 1/20], marcasLinMudas[0, 1, 1/5, 1/20]},
                    {marcasLog[0, 7],            marcasLogMudas[0, 7]}},
     GridLines -> {Join[Table[10.^k, {k, 0, 7}],
                        {{xMR, verticalEq}, {xMDE, verticalEq}}], Automatic},
     Epilog -> {
        (* labels of the equality z, next to the bottom axis (white background) *)
        Text[Framed[Style[Row[{"z\[TildeTilde]", NumberForm[zMDE, {3, 2}]}], 10, GrayLevel[0.25]],
              Background -> White, FrameStyle -> None, FrameMargins -> 1],
           Scaled[{Log[10, xMDE]/7, 0.16}], {-1.05, 0}],
        Text[Framed[Style[Row[{"z\[TildeTilde]", Round[zMR, 100]}], 10, GrayLevel[0.25]],
              Background -> White, FrameStyle -> None, FrameMargins -> 1],
           Scaled[{Log[10, xMR]/7, 0.07}], {-1.05, 0}]}],
   Placed[filaLey[{colRad, colMat, colDE},
      {"Radiaci\[OAcute]n", "Materia", "DE (CPL)"}, dash3], Below]]];
figFracCMB = figFrac[parCMB, "BAO+CMB"]
figFracCS  = figFrac[parCS,  "BAO+CMB+SNe"]


(* ==== BLOCK 9 : key numbers (check, they are printed) ==== *)
Do[Module[{p = mp[[1]], nm = mp[[2]], tH, edad, za, q0},
   tH = 9.778/p[[4]];
   edad = tH NIntegrate[1/(a EzC[1/a - 1, p[[1]], p[[2]], p[[3]], orrOf[p[[4]]]]), {a, 1*^-8, 1}];
   za = zacc[p]; q0 = qC[0, p[[1]], p[[2]], p[[3]], orrOf[p[[4]]]];
   Print[nm, ":  edad = ", NumberForm[edad, {4, 2}], " Gyr   z_acc(q=0) = ",
         NumberForm[za, {3, 2}], "   q0 = ", NumberForm[q0, {3, 2}]]],
  {mp, {{parLCDM, "LCDM"}, {parCMB, "BAO+CMB"}, {parCS, "BAO+CMB+SNe"}}}];


(* ==== BLOCK 10 -> Export of the plots to PDF. ==== *)

(* 1. We tell Mathematica to work in the same folder where this .nb file is saved. *)
SetDirectory[NotebookDirectory[]];

(* 2. We export the 7 figures to high-quality vector PDF format. *)
Export["A_Historia_Expansion_Ez.pdf", figEz];
Export["B_Desaceleracion_qz.pdf", figQz];
Export["C_Factor_Escala_at.pdf", figAt];
Export["D_Densidades_Eras_BAO_CMB.pdf", figDensCMB];
Export["D_Densidades_Eras_BAO_CMB_SNe.pdf", figDensCS];
Export["E_Fracciones_Energia_BAO_CMB.pdf", figFracCMB];
Export["E_Fracciones_Energia_BAO_CMB_SNe.pdf", figFracCS];

(* 3. Confirmation message to know that everything went well. *)
Print["\[DownExclamation]Exito! Las 7 graficas CPL-DESI se han guardado en la carpeta: ", Directory[]];


(* =====================================================================
   Notes:
   - Parameters: isolated in BLOCK 2. If your fit (ajuste_DESI_*.wl) gives you
     slightly different values for BAO+CMB or BAO+CMB+SNe, change them there.
   - The DE in CPL (red in D/E) is not flat: for these fits w was
     very negative (phantom) in the past and grows towards the future, so
     rho_DE peaks near today. Far from the data range (z>~2) it is a CPL
     extrapolation, take it with caution.
   - The dash-dot vertical lines mark the density EQUALITIES: in
     figs. D on the a axis (radiation-matter a~3e-4, matter-DE a~0.74) and in
     figs. E on the z axis (radiation-matter z~3400, matter-DE z~0.3),
     with their value next to the axis.
   - Distances: they are not compared here. The three models give D_M(z) within
     ~2-3% of each other over the data range (and the two w0waCDM fits agree
     to <0.5%), a difference invisible in a normal D(z); besides, D_A and D_L
     are only D_M rescaled by (1+z), so they add no information. What separates
     the models is seen in q(z), a(t) and the evolution of rho_DE. The distance
     plot is best left ONLY in LCDM (distances section), as an illustration
     of the definitions (D_C, D_A, D_L), the duality D_L=(1+z)^2 D_A
     and the maximum of D_A.
   ===================================================================== *)
