(* ::Package:: *)

(* =====================================================================
   plots_DE_with_w_const.wl   (Wolfram Language / Mathematica)
   Calculations and plots of the LCDM model with different energy content.

   >>> Uniform format (Frame, grid, dashed markers) set with
       SetOptions (BLOCK 3b), and labels with real symbols
       (rho, Omega, Lambda, subscripts and superscripts) on axes and legends.

   >>> Legends: they are built by hand (BLOCK 3b, function filaLey) because
       LineLegend draws the colour swatch as a glyph from an internal font
       that on some installations comes out broken (the strange symbols
       in front of each label). Drawing the swatch as a real line
       (Graphics + Line) makes the problem disappear, and putting everything
       in a single Row means the legend never splits into two lines.
   ===================================================================== *)


(* ==== BLOCK 1: constants and units ==== *)
cc  = 299792.458;        (* speed of light [km/s] *)
h   = 0.7;               (* H0 = 100 h km/s/Mpc, where I choose h=0.7, as it is a value between the two extreme ones that are responsible for the
                            Hubble tension. *)
h0  = 100 h;             (* Hubble constant today [km/s/Mpc]. *)
dH0 = cc/h0;             (* Hubble distance c/H0 [Mpc]. *)
tH  = 9.778/h;           (* Hubble time 1/H0 in [Gyr]. *)
SetDirectory[NotebookDirectory[]];


(* ==== BLOCK 2: energy content (fiducial LCDM) ==== *)
omR0 = 9.0*^-5;                  (* radiation today (photons+neutrinos), h=0.7. *)
omM0 = 0.31;                     (* matter today (baryons + dark matter). *)
omK0 = 0.0;                      (* curvature (flat universe) *)
omL0 = 1 - omM0 - omR0 - omK0;   (* Lambda from closure (~0.69). *)


(* ==== BLOCK 3: basic functions ==== *)
Ez[z_, omR_, omM_, omL_, omK_] :=
   Sqrt[omR (1 + z)^4 + omM (1 + z)^3 + omK (1 + z)^2 + omL];   (* E(z)=H/H0. *)

Ea[a_, omR_, omM_, omL_, omK_] :=
   Sqrt[omR a^-4 + omM a^-3 + omK a^-2 + omL];                  (* E as a function of a. *)

qz[z_, omR_, omM_, omL_, omK_] := Module[{s, ds},              (* q(z)=(1+z)E'/E-1. *)
   s  = omR (1 + z)^4 + omM (1 + z)^3 + omK (1 + z)^2 + omL;
   ds = 4 omR (1 + z)^3 + 3 omM (1 + z)^2 + 2 omK (1 + z);
   (1 + z) ds/(2 s) - 1];


(* ==== BLOCK 3b: common format + labels with symbols (must be evaluated). ==== *)
colUniv = {RGBColor[0.706,0.271,0.122], RGBColor[0.106,0.286,0.396],
           RGBColor[0.165,0.616,0.561], RGBColor[0.612,0.427,0.682]};
colComp = {RGBColor[0.894,0.631,0.106], RGBColor[0.106,0.286,0.396],
           RGBColor[0.706,0.271,0.122]};
(* 4 colours for the 4 distances (D_M blue, D_H purple, D_A orange, D_L green). *)
colDist = {RGBColor[0.122,0.467,0.706], RGBColor[0.580,0.404,0.741],
           RGBColor[1.000,0.498,0.055], RGBColor[0.173,0.627,0.173]};
(*  Dash patterns (solid, dashed, dotted, dash-dot) to distinguish
   the curves also in black and white / for colour-blind readers. dash3 for
   the plots with 3 curves, dash4 for those with 4. *)
dash4 = {Dashing[None], AbsoluteDashing[{9, 5}], AbsoluteDashing[{1, 4}], AbsoluteDashing[{9, 4, 1, 4}]};
dash3 = {Dashing[None], AbsoluteDashing[{9, 5}], AbsoluteDashing[{1, 4}]};
(* The styles combine colour + thickness + dash pattern. *)
estUniv = MapThread[Directive[#1, AbsoluteThickness[1.8], #2] &, {colUniv, dash4}];
estComp = MapThread[Directive[#1, AbsoluteThickness[1.8], #2] &, {colComp, dash3}];
estDist = MapThread[Directive[#1, AbsoluteThickness[1.8], #2] &, {colDist, dash4}];

(* Special lines, different from the grid. *)
lineaDisc = Directive[GrayLevel[0.4], Dashed, AbsoluteThickness[1.3]]; 
lineaPunt = Directive[GrayLevel[0.4], Dotted, AbsoluteThickness[1.5]]; 
(* dash-dot style for the equality vertical lines. *)
verticalEq = Directive[GrayLevel[0.35], AbsoluteDashing[{6, 3, 1, 3}], AbsoluteThickness[1.3]];

(* Symbol shortcuts(Greek letters).*)
sOmM = Subscript["\[CapitalOmega]", "M"]; (* Omega_M (uppercase). *)
sOmL = Subscript["\[CapitalOmega]", "\[CapitalLambda]"]; (* Omega_Lambda. *)
sOmi = Subscript["\[CapitalOmega]", "i"]; (* Omega_i. *)
sOmK = Subscript["\[CapitalOmega]", "\[Kappa]"]; (* Omega_kappa (curvature). *)
sH0  = Subscript["H", "0"]; (* H_0. *)
sRhoi= Subscript["\[Rho]", "i"]; (* rho_i. *)
sRhoc= Subscript["\[Rho]", "c,0"];(* rho_{c,0}. *)

(* Legends (lists of typeset labels).*)
legUniv = {Row[{"EdS (", sOmM, "=1)"}], "\[CapitalLambda]CDM",
           Row[{"Abierto (", sOmM, "=0.3, ", sOmK, "=0.7)"}], Row[{"De Sitter (", sOmL, "=1)"}]};
legComp2 = {Row[{"Radiaci\[OAcute]n \[Proportional] ", Superscript["a", "-4"]}],
            Row[{"Materia \[Proportional] ", Superscript["a", "-3"]}],
            "\[CapitalLambda] = cte"};
legComp3 = {"Radiaci\[OAcute]n", "Materia", "\[CapitalLambda]"};
legDist  = {Row[{Subscript["D", "C"], " comovil"}],
            Row[{Subscript["D", "A"], " diam. angular"}],
            Row[{Subscript["D", "L"], " luminosidad"}]};
(* Short legend with the 4 distances (same order as the Plot of fig.6). *)
legDistCorto = {Subscript["D", "M"], Subscript["D", "H"], Subscript["D", "A"], Subscript["D", "L"]};

(* -- Manual legend: the swatch is a real line (not a glyph), so it
      cannot come out broken; all of it goes in one Row, so it never splits. -- *)
(* muestra accepts a DASH pattern (dash) besides the colour, so the
   legend shows the same line pattern as the curve (distinguishable in B&W). *)
muestra[col_, dash_ : Dashing[None]] := Graphics[
   {col, AbsoluteThickness[2.4], dash, Line[{{0, 0}, {1, 0}}]},
   ImageSize -> {30, 10}, AspectRatio -> Full,
   PlotRange -> {{0, 1}, {-1, 1}}, PlotRangePadding -> None,
   ImagePadding -> None, BaselinePosition -> Center,
   Frame -> False, Axes -> False, GridLines -> None,   (* <- no frame and no grid *)
   Prolog -> {}, Epilog -> {}, Background -> None];

itemLey[col_, lab_, dash_ : Dashing[None]] := Row[{muestra[col, dash], Spacer[5], Style[lab, Black, 12]}];

filaLey[cols_, labs_] := Row[Riffle[MapThread[itemLey, {cols, labs}], Spacer[22]]];
(* Variant with a list of dash patterns (one per curve, in the same order). *)
filaLey[cols_, labs_, dashes_] := Row[Riffle[MapThread[itemLey, {cols, labs, dashes}], Spacer[22]]];

(* Ticks on the logaritmic axes, all as 10^k ----------------------
   By default Mathematica mixes on the same logarithmic axis the decimal
   notation (0.001, 0.1, 10, 1000) with the scientific one (10^-5, 10^5, ...).
   Here the ticks are built BY HAND so that all of them come out as 10^k.

     marcasLog[kmin, kmax]        -> label 10^k at every decade
     marcasLog[kmin, kmax, paso]  -> label only every 'paso' decades; the
                                     ones in between stay as short ticks.
                                     withou label (so text does not pile up).
     marcasLogMudas[...]          -> the same ticks without labels, for the
                                     top and right sides of the frame. *)
marcasLog[kmin_Integer, kmax_Integer, paso_Integer : 1] :=
  Module[{conEtiq, sinEtiq},
   conEtiq = Range[kmin, kmax, paso];                  (* labelled decades *)
   sinEtiq = Complement[Range[kmin, kmax], conEtiq];   (* intermediate decades *)
   Join[
     Table[{10.^k, Superscript[10, k], {0.010, 0}}, {k, conEtiq}],
     Table[{10.^k, "", {0.005, 0}}, {k, sinEtiq}]]];

marcasLogMudas[kmin_Integer, kmax_Integer, paso_Integer : 1] :=
  ({#[[1]], "", #[[3]]} &) /@ marcasLog[kmin, kmax, paso];

(* -- variant for log axes with few decades: label 10^k
      at every decade and also short ticks without label at 2,3,...,9 of
      each decade, which is the usual look of a logarithmic axis. Ticks
      falling outside the PlotRange are simply not drawn.                      *)
marcasLogFinas[kmin_Integer, kmax_Integer] :=
  Join[
    Table[{10.^k, Superscript[10, k], {0.010, 0}}, {k, kmin, kmax}],
    Flatten[Table[{m 10.^k, "", {0.005, 0}}, {k, kmin - 1, kmax}, {m, 2, 9}], 1]];

marcasLogFinasMudas[kmin_Integer, kmax_Integer] :=
  ({#[[1]], "", #[[3]]} &) /@ marcasLogFinas[kmin, kmax];

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

(* -- common default format for ALL the plots -- *)
SetOptions[{Plot, LogPlot, LogLogPlot, LogLinearPlot, ParametricPlot},
   (* LogPlot is added to the list, because fig.6 (distances) has
        the y axis on a logarithmic scale and would otherwise be left without
        this common format (frame, grid, typography, size...).                 *)
   Frame -> True, Axes -> False, GridLines -> Automatic,
   GridLinesStyle -> Directive[GrayLevel[0.85], AbsoluteThickness[0.4]],
   FrameStyle -> Directive[Black, AbsoluteThickness[1]],
   LabelStyle -> Directive[Black, 12], Background -> White,
   ImageSize -> 480, AspectRatio -> 0.72, PlotRangePadding -> Scaled[0.015],
   ImagePadding -> {{Automatic, 22}, {Automatic, Automatic}}
   (* Fixed right margin (22 pt) so that the right line of the frame is not
        cut off in LogLinearPlot/ParametricPlot; the rest, automatic. *)
];


(* ==== BLOCK 4: key numbers (verification) ==== *)
zEqMR = omM0/omR0 - 1;
zEqML = (omL0/omM0)^(1/3) - 1;
zAcc  = z /. FindRoot[qz[z, omR0, omM0, omL0, omK0] == 0, {z, 0.6}];
q0    = qz[0, omR0, omM0, omL0, omK0];
age   = tH NIntegrate[1/(a Ea[a, omR0, omM0, omL0, omK0]), {a, 1*^-8, 1}];
Print["z_eq(mat-rad)     = ", zEqMR];
Print["z_eq(mat-Lambda)  = ", zEqML];
Print["z_acc (q=0)       = ", zAcc];
Print["q0                = ", q0, "   (Or+Om/2-OL = ", omR0 + omM0/2 - omL0, ")"];
Print["Edad del universo = ", age, " Gyr"];


(* ==== BLOCK 5 -> FIGURE 1: E(z) ==== *)
fig1 = Legended[
  Plot[{
     Ez[z, 0,    1.0,  0.0, 0.0],
     Ez[z, omR0, omM0, omL0, 0.0],
     Ez[z, 0,    0.3,  0.0, 0.7],
     Ez[z, 0,    0.0,  1.0, 0.0]
    }, {z, 0, 3},
    PlotStyle  -> estUniv,
    FrameLabel -> {"z", Row[{"E(z) = H/", sH0}]},
    PlotRange  -> All],
  Placed[filaLey[colUniv, legUniv, dash4], Below]]


(* ==== BLOCK 6 -> FIGURE 2: densities and eras rho_i(z)/rho_c,0(LCDM) ==== *)
(* The two equality vertical lines are in dash-dot style
   (verticalEq, as in graficas_CPL_DESI) and with a LABEL of the value of a.
   Positions (a axis):
     radiation-matter:   Or a^-4 = Om a^-3  ->  a = omR0/omM0.
     matter-Lambda:      Om a^-3 = omL0     ->  a = (omM0/omL0)^(1/3).          *)
fig2 = Legended[
  LogLogPlot[{omR0 a^-4, omM0 a^-3, omL0}, {a, 1*^-6, 100},
    PlotStyle  -> estComp,
    FrameLabel -> {"a", Row[{sRhoi, " / ", sRhoc}]},
    PlotRange  -> {{1*^-6, 100}, {1*^-12, 1*^20}},
    (* Ticks of both axes all as 10^k (without mixing 0.001 with 10^-5).
       y axis: 10^-12 ... 10^20, with a label every 3 decades (like the grid).
       x axis: 10^-6 ... 10^2, with a label at every decade.                    *)
    FrameTicks -> {{marcasLog[-12, 20, 3], marcasLogMudas[-12, 20, 3]},
                   {marcasLog[-6,  2],    marcasLogMudas[-6,  2]}},
    GridLines  -> {
        Join[ Table[10.^k, {k, -6, 2}],
              {{omR0/omM0, verticalEq}, {(omM0/omL0)^(1/3), verticalEq}} ],
        Table[10.^k, {k, -12, 18, 3}]
    },
    Epilog -> {   (* Labels of the equality a (white background); x fraction of the log axis [1e-6,100] = (Log10[a]+6)/8 *)
       Text[Framed[Style[Row[{"a\[TildeTilde]", ScientificForm[omR0/omM0, 2]}], 10, GrayLevel[0.25]],
             Background -> White, FrameStyle -> None, FrameMargins -> 1],
          Scaled[{(Log[10, omR0/omM0] + 6)/8, 0.06}], {-1.2, 0}],
       Text[Framed[Style[Row[{"a\[TildeTilde]", NumberForm[(omM0/omL0)^(1/3), {3, 2}]}], 10, GrayLevel[0.25]],
             Background -> White, FrameStyle -> None, FrameMargins -> 1],
          Scaled[{(Log[10, (omM0/omL0)^(1/3)] + 6)/8, 0.06}], {1.2, 0}]}],
  Placed[filaLey[colComp, legComp2, dash3], Below]]


(* ==== BLOCK 7 -> FIGURE 3: Omega_i(z)=rho_i(z)/rho_c(z) (LCDM) ==== *)
(* Equality vertical lines in verticalEq + label of the value of z.
   Positions (1+z axis):
     radiation-matter:   1+z = omM0/omR0.
     matter-Lambda:      1+z = (omL0/omM0)^(1/3).                               *)
fig3 = Legended[
  LogLinearPlot[
    Evaluate@With[{s = omR0 u^4 + omM0 u^3 + omL0},
       {omR0 u^4/s, omM0 u^3/s, omL0/s}],
    {u, 1, 1*^7},
    PlotStyle  -> estComp,
    PlotRange  -> {{1, 1*^7}, {-0.02, 1.02}},
    (* X axis (log) all as 10^k: 10^0, 10^1, ..., 10^7, instead of
       1, 10, 100, 1000, 10^4, ... (which mixes decimal and scientific).
       The y axis is linear and stays as always: 0.0, 0.2, ..., 1.0.           *)
    FrameTicks -> {{marcasLin[0, 1, 1/5, 1/20], marcasLinMudas[0, 1, 1/5, 1/20]},
                   {marcasLog[0, 7],            marcasLogMudas[0, 7]}},
    GridLines  -> {                                   (* As in fig.2: grid + eras. *)
       Join[ Table[10.^k, {k, 0, 7}],                  (* Vertical grid (decades). *)
             {{omM0/omR0, verticalEq},                 (* Radiation-matter intersection (dash-dot). *)
              {(omL0/omM0)^(1/3), verticalEq}} ],       (* Matter-Lambda intersection (dash-dot). *)
       Automatic                                       (* Automatic horizontal grid. *)
    },
    FrameLabel -> {"1+z", Row[{sOmi, "(z)"}]},
    Epilog -> {   (* Labels of the equality z (white background); x fraction of the log axis [1,1e7] = Log10[1+z]/7. *)
       Text[Framed[Style[Row[{"z\[TildeTilde]", NumberForm[(omL0/omM0)^(1/3) - 1, {3, 2}]}], 10, GrayLevel[0.25]],
             Background -> White, FrameStyle -> None, FrameMargins -> 1],
          Scaled[{Log[10, (omL0/omM0)^(1/3)]/7, 0.16}], {-1.05, 0}],
       Text[Framed[Style[Row[{"z\[TildeTilde]", Round[omM0/omR0 - 1, 100]}], 10, GrayLevel[0.25]],
             Background -> White, FrameStyle -> None, FrameMargins -> 1],
          Scaled[{Log[10, omM0/omR0]/7, 0.07}], {-1.05, 0}]}],
  Placed[filaLey[colComp, legComp3, dash3], Below]]


(* ==== BLOCK 8 -> FIGURE 4: a(t) ==== *)
tofa[a_?NumericQ, omR_, omM_, omL_, omK_] :=
   tH NIntegrate[1/(x Ea[x, omR, omM, omL, omK]), {x, 1, a}];

fig4 = Legended[
  ParametricPlot[{
     {tofa[a, 0,    1.0,  0.0, 0.0], a},
     {tofa[a, omR0, omM0, omL0, 0.0], a},
     {tofa[a, 0,    0.3,  0.0, 0.7], a},
     {tofa[a, 0,    0.0,  1.0, 0.0], a}
    }, {a, 0.001, 4.4},   (* a from 0.001 to 4.4: at the bottom the curves with matter go down to the axis (a->0, Big Bang)
     and at the top the fastest one (de Sitter, a~4.2) reaches t=20, without being cut off at the top or at the bottom. *)
    PlotPoints -> 120, MaxRecursion -> 4,   (* Fine sampling: draws the whole almost vertical stretch of the Big Bang until it touches a=0 
    (otherwise ParametricPlot left it short at the bottom). *)
    PlotStyle  -> estUniv,
    FrameLabel -> {Row[{"t - ", Subscript["t", "0"], " [Gyr]  (0 = hoy)"}], "a(t)"},
    PlotRange  -> {{-15, 20}, {0, 4.3}},   (* Y axis up to 4.3 to see all the curves up to t=20 without them being cut off at the top. *)
    Epilog     -> {
       {lineaPunt, Line[{{-15, 1}, {20, 1}}]},    (* a = 1. *)
       {lineaPunt, Line[{{0, 0}, {0, 4.3}}]}        (* t - t0 = 0 *) (* Vertical line up to the new y maximum. *)
    }],
  Placed[filaLey[colUniv, legUniv, dash4], Below]]


(* ==== BLOCK 9 -> FIGURE 5: q(z) ==== *)
(* Dash-dot vertical line (verticalEq) at z_acc (where LCDM crosses q=0,
   i.e. the onset of acceleration) + label with that z. The other three
   models do not cross q=0 (EdS q=1/2, de Sitter q=-1, Open q>0), so the
   vertical line marks the z_acc of LCDM.  I fix PlotRange so that the
   vertical line is not cut off.                                               *)
fig5 = Legended[
  Plot[{
     qz[z, 0,    1.0,  0.0, 0.0],
     qz[z, omR0, omM0, omL0, 0.0],
     qz[z, 0,    0.3,  0.0, 0.7],
     qz[z, 0,    0.0,  1.0, 0.0]
    }, {z, 0, 3},
    PlotStyle  -> estUniv,
    FrameLabel -> {"z", "q(z)"},
    PlotRange  -> {{0, 3}, {-1.05, 0.55}},
    Epilog     -> {
       {lineaDisc, Line[{{0, 0}, {3, 0}}]},                (* q=0. *)
       {verticalEq, Line[{{zAcc, -1.05}, {zAcc, 0.55}}]},  (* Vertical line at z_acc (LCDM). *)
       {Black, PointSize[0.018], Point[{zAcc, 0}]},        (* Onset of acceleration. *)
       Text[Framed[Style[Row[{"z\[TildeTilde]", NumberForm[zAcc, {3, 2}]}], 10, GrayLevel[0.25]],
             Background -> White, FrameStyle -> None, FrameMargins -> 1],
          {zAcc, -0.9}, {-1.05, 0}]                        (* Label of z_acc *)
    }],
  Placed[filaLey[colUniv, legUniv, dash4], Below]]


(* ==== BLOCK 10 -> FIGURE 6: distances (LCDM) ==== *)
dCf[z_?NumericQ] := dH0 NIntegrate[1/Ez[zp, omR0, omM0, omL0, omK0], {zp, 0, z}];

(* For D_H = c/H(z) = dH0/E(z) (Hubble distance, the RADIAL BAO
   observable), it is not an integral like the other three; it
   is local (D_H = dD_C/dz) and decreases with z (H grows with z). Order/colours
   of fig.6: D_M blue, D_H purple, D_A orange, D_L green. Recall (flat universe):
   D_A = D_M/(1+z),  D_L = (1+z) D_M  (so D_A and D_L come from D_M).           *)
(* Y axis on a logarithmic scale (LogPlot instead of Plot).
   Reason: at z=3 the four distances are D_M=6293, D_H=945, D_A=1573 and
   D_L=25172 Mpc. With a linear y axis, the automatic PlotRange discarded the
   upper part of D_L (it was cut off around z~1.8); and forcing PlotRange->All
   so that it fitted, the maximum of D_A (1734 Mpc) stayed in the bottom 7% of
   the height, i.e. the three small distances were squashed against the axis.
   In log scale the data range (42 to 25172 Mpc) is only ~2.8 decades:
   all four curves fit entirely and remain well separated from each other,
   including the maximum of D_A and the crossing D_H = D_A.                    *)
fig6 = Legended[
  LogPlot[{dCf[z], dH0/Ez[z, omR0, omM0, omL0, omK0], dCf[z]/(1 + z), (1 + z) dCf[z]},
    {z, 0.01, 3},
    PlotStyle  -> estDist,
    FrameLabel -> {"z", "distancia [Mpc]"},
    PlotRange  -> {{0, 3}, {30., 30000.}},   (* exactly 3 decades: 3*10^1 ... 3*10^4. *)
    (* ticks: y axis as 10^k (10^2, 10^3, 10^4) with short ticks at 2..9 of
       each decade; linear x axis as before (0.0, 0.5, ..., 3.0).              *)
    FrameTicks -> {{marcasLogFinas[2, 4],         marcasLogFinasMudas[2, 4]},
                   {marcasLin[0, 3, 1/2, 1/10],   marcasLinMudas[0, 3, 1/2, 1/10]}},
    (* grid set by hand so that it falls exactly where the labels are:
       vertical lines every 0.5 in z, horizontal ones at the decades (and the
       intermediate lines 2..9 of each decade, lighter, to read the log axis). *)
    GridLines  -> {Range[0, 3, 1/2],
                   Join[Table[10.^k, {k, 2, 4}],
                        Flatten[Table[{m 10.^k,
                              Directive[GrayLevel[0.93], AbsoluteThickness[0.4]]},
                            {k, 1, 4}, {m, 2, 9}], 1]]}],
  Placed[filaLey[colDist, legDistCorto, dash4], Below]]

(* =====================================================================
   Notes:
   - The format of all Six is controlled in the SetOptions of BLOCK 3b.
   - The Legend is drawn by hand with filaLey (BLOCK 3b): the colour
     swatch is a real line (Graphics+Line), not a glyph, and everything goes
     in one Row, which is why it does not break or split into two lines. The
     swatch size is adjusted in 'muestra' (ImageSize) and the spacing between
     entries in 'filaLey' (Spacer[22]).
   - Greek letters: "\[Rho]"=rho, "\[CapitalOmega]"=Omega, "\[CapitalLambda]"=
     Lambda, "\[Proportional]"= the proportionality symbol. On evaluation,
     Mathematica turns them into the symbol (you will see real rho, Omega...).
   - To vary the energy content change omM0/omL0/omK0 (BLOCK 2).
   - The equality vertical lines (fig.2, fig.3) and the q=0 one (fig.5)
     are in dash-dot style 'verticalEq' with their value on the axis (a in
     fig.2, z in fig.3 and fig.5). Those of fig.2/3 are the SAME intersections,
     one on the a axis and the other on the 1+z axis (a = 1/(1+z)).
   - Fig.6 now has 4 distances, including D_H = c/H(z) (Hubble /
     radial BAO). D_H DECREASES with z (it starts at c/H0 ~ 4283 Mpc) while the
     others grow; the maximum of D_A falls exactly where D_H = D_A (~z=1.6).
   - Fig.6 uses a logarithmic y axis (LogPlot). D_L grows up to
     25172 Mpc at z=3 while D_A stays at ~1700: on a linear scale either D_L
     was cut off or the other three were squashed. If the linear scale were
     still preferred, just change LogPlot to Plot, remove the FrameTicks and
     the GridLines of that figure and set PlotRange -> {{0, 3}, {0, 26000}}
     (but then D_H, D_M and D_A end up squeezed against the axis).
   ===================================================================== *)


(* ==== BLOCK 11 -> EXPORT OF THE PLOTS TO PDF ==== *)

(* 1. We tell Mathematica to work in the same folder where this .nb file is saved. *)
SetDirectory[NotebookDirectory[]];

(* 2. We export the 6 variables to high-quality vector PDF format. *)
Export["1_Historia_Expansion_Ez.pdf", fig1];
Export["2_Densidades_Eras.pdf", fig2];
Export["3_Fracciones_Energia_Omega.pdf", fig3];
Export["4_Factor_Escala_at.pdf", fig4];
Export["5_Desaceleracion_qz.pdf", fig5];
Export["6_Distancias_Cosmologicas.pdf", fig6];

(* 3. Confirmation message to know that everything went well. *)
Print["\[DownExclamation]Exito! Las 6 graficas se han guardado en la carpeta: ", Directory[]];
