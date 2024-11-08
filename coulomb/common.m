(* the operator referred to as \mathcal{D} in the text *)

ζ[r_] = (I*(r*D[#,{r,1}] + 5/2*#))&;

(* \mathcal{D}^*, the adjoint operator of \mathcal{D} *)

ζ2[r_] = (I*(r*D[#,{r,1}] + 1/2*#))&;

(* apply a polynomial p(op) to a function f *)

PolyOp[p_,op_,f_] := Sum[Coefficient[p[s], s, j]*Nest[op,f,j], {j,0,Exponent[p[s],s]}];

(* Meixner-Pollaczek polynomial *)

Pmp[n_, x_, λ_, φ_] := Pochhammer[2*λ,n]/n! * Exp[I*n*φ] * Hypergeometric2F1[-n, λ + I*x, 2*λ, 1 - Exp[-2*I*φ]];

hmp[n_, λ_, φ_] := 2*π*Gamma[n+2*λ]/((2*Sin[φ])^(2*λ) * n!); (* normalisation constant *)

kmp[n_, λ_, φ_] := (2*Sin[φ])^n/(n!); (* coeff of x^n *)

(* symmetric continuous Hahn polynomial *)

pn[n_,a_,b_,x_] := I^n * Pochhammer[2*a,n]*Pochhammer[a+b,n]/(n!) * HypergeometricPFQ[{-n,n+2a+2b-1,a+I*x},{2a,a+b},1];

hn[n_,a_,b_] := 2*Pi*Gamma[n+2*a]*Gamma[n+2*b]*(Gamma[n+a+b])^2/((2*n+2*a+2*b-1)*Gamma[n+2*a+2*b-1]*(n!)); (* normalisation const for continuous Hahn polys *)

kn[n_,a_,b_] := Pochhammer[n + 2*a + 2*b - 1,n]/(n!);

βn[n_,a_,b_] := n*(n+2a-1)*(n+2b-1)*(n+2a+2b-2)/(4*(2n+2a+2b-3)*(2n+2a+2b-1));

(* Laguerre polynomial with non-negative degree *)

LagL[n_,a_,x_] := If[n < 0, 0, LaguerreL[n,a,x]];

(* reverse Bessel polynomial *)

θn[n_,x_] := (n!)/(-2)^n*LagL[n,-2*n-1,2*x];


