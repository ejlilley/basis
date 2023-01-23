(* Mathematica code to produce symbolic expressions/plots for the
isochrone-adapted basis set *)

(* Note: This code uses the Modified Chebyshev algorithm (see Sec. 5.2
in the paper), whereas the julia code (in Isochrone.jl) uses the
Discretized Stieltjes procedure (see Sec. 5.1). The code here mostly
serves to check that the two algorithms give the same result
(numerical instability notwithstanding). *)

Dop[r_] = (I*(r*D[#,{r,1}] + 5/2*#))&; (* the operator referred to in the paper as "D" (Eq. (2.4)) *)

Dop2[r_] = (I*(r*D[#,{r,1}] + 1/2*#))&; (* D - 2i *)

(* The l-dependent potential, density and polynomial weight function. See Eq. (5.1). *)

phi0l[l_,r_] := r^l/(1 + Sqrt[1 + r^2])^(1+2*l);

rho0l[l_, r_] := FullSimplify[(x^(-2)*D[x^2*D[phi0l[l,x],x],x] - l*(l+1)/x^2*phi0l[l,x]), Assumptions -> {x≥0,Element[x,Reals],l≥0,Element[l,Integers]}] /. x -> r;

ωl[l_, s_] := Re[((l+1/2)^2 + s^2)/(8*Pi^2) * MellinTransform[phi0l[l,r], r, 1/2 + I*s] * MellinTransform[phi0l[l,r], r, 1/2 - I*s]];

(* Exact expressions for the isochrone modified moments. See Appendix G. *)

fjl[j_,l_,y_] := 2^(l+3/2)*(3l+5/2)*Pochhammer[3l+7/2,j]/(j! * Pochhammer[2l+2,j])*(1-x)^j*Nest[(1-x)*D[(1-x)^(-1)*#,x]&,(1-x)^(l+3/2)/(x^(3l+5/2))*Beta[x,3l+5/2,-l-1/2],j] /. x -> y;

Clear[Kjl]; Kjl[j_,l_] := Kjl[j,l] = Simplify[(1 + 2l)*Gamma[l+3/2]*Gamma[2l+j+2]/(2^(3l+5/2)*Gamma[3l+j+7/2])*fjl[j,l,1/2], Assumptions -> {l>=0,Element[l,Integers]}];

Clear[qjkl]; qjkl[0,0,l_] := 1; qjkl[1,0,l_] := -I/2*(1+2l); qjkl[1,1,l_] := I*(1+2*l); qjkl[j_,k_,l_] := qjkl[j,k,l] = FullSimplify[If[k<0,0,If[k>2j-1,0,I*((1+2l)*qjkl[j-1,k-1,l] - (1/2+l+k)*qjkl[j-1,k,l] + (k-2)*qjkl[j-1,k-2,l]) - βtkl[j-1,l]*qjkl[j-2,k,l]]], Assumptions -> {l>=0,Element[l,Integers],j>=0,Element[j,Integers]}];

 (* The (monic) Hermite polynomials *)

He[n_, x_] := 2^(-n/2)*HermiteH[n,x/Sqrt[2]];

ptjl[j_, l_, s_] := He[j,s];

βtkl[k_, l_] := k; (* the recurrence coeff for the monic Hermite polynomials Eq (5.9) *)

μtjl[j_, l_] := If[OddQ[j], 0, 1/(4*Pi)*Sum[qjkl[j,k,l]*(Kjl[k+1,l]+2*(l+1)*Kjl[k,l]),{k,0,Max[2*j-1,0]}]]; (* the modified moments Eq (5.8) *)

Clear[σjkl]; σjkl[-1, k_, l_] := 0; σjkl[0, k_, l_] := μtjl[k, l]; σjkl[j_, k_, l_] := σjkl[j, k, l] = FullSimplify[σjkl[j-1,k+1,l] - βjl[j-1,l]*σjkl[j-2,k,l] + βtkl[k,l]*σjkl[j-1,k-1,l]]; (* the mixed moments Eq ((5.10) *)

Nnl[n_,l_] := 4*Pi*σjkl[n,n,l]; (* The normalisation constant *)

Clear[βjl]; βjl[-1, l_] := 0; βjl[0, l_] := FullSimplify[μtjl[0, l]]; βjl[j_, l_] := FullSimplify[σjkl[j,j,l]/σjkl[j-1,j-1,l]]; (* The recurrence coeff Eq (5.11) *)

Clear[pnl]; pnl[-1, l_, s_] := 0; pnl[0, l_, s_] := 1; pnl[n_, l_, s_] := pnl[n, l, s] = Expand[Expand[s * pnl[n - 1, l, s]] - Expand[βjl[n-1,l] pnl[n - 2, l, s]]]; (* The index-raising polynomials *)

ρnl[n_, l_, r_] := Sum[Coefficient[I^(-n)*pnl[n,l,s], s, j]*Nest[Dop[r],rho0l[l,r],j], {j,0,n}];

Φnl[n_, l_, r_] := Sum[Coefficient[I^(-n)*pnl[n,l,s], s, j]*Nest[Dop2[r],phi0l[l,r],j], {j,0,n}];

(* The basis functions (ρnl, Φnl) and the recurrence coeff βjl can be
evaluated symbolically or numerically, which is mostly useful to check
that the julia code (in Isochrone.jl) is correct. It is important to
tell Mathematica to use enough digits of accuracy internally when
evaluating the (extremely complicated) symbolic expressions. For
example, to evaluate the recurrence coeff for (n=10,l=5), the
following is needed:

   N[Simplify[βjl[10,5]], 50]

   *)


(* These commands produce the CSV tables used in the plots of the
isochrone basis functions in the paper (Fig. 1) *)

phitable = Table[FullSimplify[Φnl[n,0,r]/Sqrt[Nnl[n,0]], Assumptions -> {Element[r,Reals],r>0}], {n,0,3}];

phitable2 = Table[FullSimplify[Φnl[n,1,r]/Sqrt[Nnl[n,1]], Assumptions -> {Element[r,Reals],r>0}], {n,0,3}];

d = Table[N[{r, phitable[[1]],phitable[[2]],phitable[[3]],phitable[[4]]} /. r -> 10^t, 50], {t,-3,3,0.05}]; Export["isochrone_potentials.dat", d, "CSV"];

d2 = Table[N[{r, phitable2[[1]],phitable2[[2]],phitable2[[3]],phitable2[[4]]} /. r -> 10^t, 50], {t,-3,3,0.05}]; Export["isochrone_potentials_l1.dat", d2, "CSV"];


