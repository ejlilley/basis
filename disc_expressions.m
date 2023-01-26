(* Mathematica code to produce symbolic expressions for the thin-disc
basis sets mentioned in the paper *)

ξ[R_] = (I*(R*D[#,{R,1}] + 3/2*#))&; (* The operator referred to in the paper as "A" (see Eq. (2.5)) *)

ξ2[R_] = (I*(R*D[#,{R,1}] + 1/2*#))&; (* the previous operator, minus i *)

Km[m_, s_] := (Gamma[(m + 3/2 + s)/2]*Gamma[(m + 3/2 - s)/2])/(Gamma[(m + 1/2 + s)/2]*Gamma[(m + 1/2 - s)/2]); (* K_m is defined in Eq. (C.7) *)

(* Meixner-Pollaczek polynomials (see App. E.1) *)

Pmp[n_, x_, λ_, φ_] := Pochhammer[2*λ,n]/n! * Exp[I*n*φ] * Hypergeometric2F1[-n, λ + I*x, 2*λ, 1 - Exp[-2*I*φ]];

kmp[n_, λ_, φ_] := 2*π*Gamma[n+2*λ]/((2*Sin[φ])^(2*λ) * n!); (* normalisation constant *)

(* Symmetric continuous Hahn polynomials (see App. E.2) *)

pnab[n_,a_,b_,x_] := I^n*Pochhammer[2*a,n]*Pochhammer[a+b,n]/(n!)*HypergeometricPFQ[{-n,n+2a+2b-1,a+I*x/2},{2*a,a+b},1]; (* note argument x is halved *)

hn[n_,a_,b_] := 2*Pi*Gamma[n+2*a]*Gamma[n+2*b]*(Gamma[n+a+b])^2/((2*n+2*a+2*b-1)*Gamma[n+2*a+2*b-1]*(n!)); (* normalisation constant *)

(***********************************************************)

(* Clutton-Brock's 1972 basis set, using the Aoki & Iye normalisation (see Sec. 4.2.1) *)

psinm[n_, m_, R_] := -R^m/(1 + R^2)^(m+1/2) * GegenbauerC[n,m+1/2,(R^2 - 1)/(R^2 + 1)]; (* Eq. (4.12) *)

sigmanm[n_, m_, R_] := (n + m + 1/2)/(2*π) * (-psinm[n,m,R]/(1 + R^2));

Nnm[n_,m_] := -π*Gamma[n+2m+1]/(2^(4m+2) * n! * Gamma[m+1/2]^2);

Pnm[n_,m_,s_] := I^n*Pmp[n,s,m+1/2,π/2]; (* index-raising polynomial for Clutton-Brock's basis set (Eq. (4.11)) *)

Ωm[m_,s_] := FullSimplify[MellinTransform[sigmanm[0,m,R],R,3/2+I*s]*MellinTransform[sigmanm[0,m,R],R,3/2-I*s]/(4*π*Km[m, I*s])]; (* the polynomial weight function Eq. (4.14) *)

Pξsigma[n_,m_,R_] := Sum[Coefficient[Pnm[n,m,s],s,j]*Nest[ξ[x],sigmanm[0,m,x],j], {j,0,n}] /. x -> R; (* expressing the density basis functions via their index-raising-polynomial representation *)

(* Now check the correctness of the result by evaluating e.g.
   FullSimplify[Pξsigma[n,m,R]/sigmanm[n,m,R]]
   which should equal 1 for all values of n,m *)

(***********************************************************)

(* Qian's general-k family, using the same normalisation as Qian (Sec. 4.2.2) *)

sigmaknm[k_,-1,m_,R_] := 0; sigmaknm[k_,0,m_,R_] := Sqrt[Pi]*Gamma[m+k+3/2]/Gamma[m+k+1] * R^m/(1 + R^2)^(m+k+3/2); sigmaknm[0,n_,m_,R_] := Sqrt[Pi]*Gamma[n+m+3/2]/Gamma[n+2m+1]*Pochhammer[1/2,m]*4^m*R^m/(1 + R^2)^(3/2+m)*GegenbauerC[n,m+1/2,(R^2-1)/(R^2+1)]; sigmaknm[k_,n_,m_,R_] := sigmaknm[k,n,m,R] = (n+m-1/2)/(n+2m+2k)*sigmaknm[k,n-1,m,R] + (n+m+2k-1/2)/(n+2m+2k)*sigmaknm[k-1,n,m,R] - (n+1)/(n+2m+2k)*sigmaknm[k-1,n+1,m,R]; (* Qian's explicit(!) expressions for the density basis functions *)

psinmk[k_,-1,m_,R_] := 0; psinmk[k_,0,m_,R_] := π * Beta[m+1/2,1/2] * R^m * (1 + R^2)^(-m-1/2) * Hypergeometric2F1[-k,m+1/2,m+1,R^2/(1 + R^2)]; psinmk[0,n_,m_,R_] := Pi^(3/2)*Gamma[n+m+1/2]/(m! * Pochhammer[2m+1,n]) * R^m/(1 + R^2)^(m+1/2) * GegenbauerC[n,m+1/2,(R^2-1)/(R^2+1)]; psinmk[k_,n_,m_,R_] := psinmk[k,n,m,R] = (n+m-1/2)/(n+2m+2k)*psinmk[k,n-1,m,R] + (n+m+2k-1/2)/(n+2m+2k)*psinmk[k-1,n,m,R] - (n+1)/(n+2m+2k)*psinmk[k-1,n+1,m,R]; (* and same for the potential (but rewriting the zeroth order to use a 2F1 function) *)

Nknm[k_,n_,m_] := Pi/(4 * (Gamma[1 + k + m])^2) * hn[n,1/4+m/2,3/4+m/2+k]/(Pochhammer[k+m+1,n])^2; (* normalisation constant (Eq. (4.18)) *)

Pnmk[k_,n_,m_,s_] := I^(n)/(Pochhammer[k+m+1,n]) * pnab[n,1/4+m/2,3/4+m/2+k,s]; (* index-raising polynomial for Qian's basis set (Eq. (4.17)) *)

Ωmk[k_,m_,s_] := FullSimplify[MellinTransform[sigmaknm[k,0,m,R],R,3/2+I*s]*MellinTransform[sigmaknm[k,0,m,R],R,3/2-I*s]/(4*π*Km[m, I*s])]; (* the polynomial weight function *)

Pξsigmaknm[k_,n_,m_,R_] := Sum[Coefficient[Pnmk[k,n,m,s], s, j]*Nest[ξ[r],sigmaknm[k,0,m,r],j], {j,0,n}] /. r -> R; (* expressing density basis function via the index-raising polynomial *)

(* Now can check that the result is correct by evaluating e.g. 

   FullSimplify[Pξsigmaknm[k,n,m,R]/sigmaknm[k,n,m,R]]

   for various values of (k,n,m), and checking that the result is 1 *)


(***********************************************************)

(* Not including code for Qian's Gaussian family as there appear to be errors in his original paper. *)


(***********************************************************)

(* Exponential disc (Sec. 4.2.4) *)

sigma0me[m_, R_] := R^m * Exp[-R]; (* exponential disc density with correct small-R multipole dependence for m>0 (Eq. (4.22)) *)

psi0metest[m_,R_] := -(2^(1 + m)*Sqrt[Pi]*MeijerG[{{(1 - m)/2}, {}}, {{m/2, (2 + m)/2}, {-1/2*m}}, R^2/4]); (* The potential as spat out by Mathematica when performing the requisite Hankel transforms *)

psi0me[0,R_] := -π*R*(BesselI[0,R/2]*BesselK[1,R/2] - BesselI[1,R/2]*BesselK[0,R/2]); psi0me[m_,R_] := psi0me[m,R] = FullSimplify[x*D[psi0me[m-1,x], {x,2}] + (1-2*m)*D[psi0me[m-1,x],x] + (m^2-1)/x*psi0me[m-1,x]] /. x -> R; (* the potential as calculated using the result derived in the paper (Eq. (4.25)) *)

Ωme[m_,s_] := FullSimplify[MellinTransform[sigma0me[m,R],R,3/2+I*s]*MellinTransform[sigma0me[m,R],R,3/2-I*s]/(4*π*Km[m, I*s])]; (* the polynomial weight function *)

Pnme[n_, m_, s_] := pnab[n,m/2+1/4,m/2+5/4,s]; (* index-raising polynomial (Eq. (4.23)) *)

Nnme[n_, m_] := 2*Pi*2^(2*m)/Pi^2*hn[n,m/2+1/4,m/2+5/4]; (* normalisation constant *)

(* Now express the higher-order basis functions using the index-raising polynomial representation *)

Pnmeξ[n_,m_,R_,f_] := Sum[Coefficient[Pnme[n,m,s],s,j]*Nest[ξ[R],f,j], {j,0,n}];

Pnmeξ2[n_,m_,R_,f_] := Sum[Coefficient[Pnme[n,m,s],s,j]*Nest[ξ2[R],f,j], {j,0,n}];

sigmanme[n_, m_, R_] := I^(-n)*FullSimplify[Pnmeξ[n,m,x,sigma0me[m,x]]] /. x -> R;

psinme[n_, m_, R_] := I^(-n)*Simplify[Pnmeξ2[n,m,x,psi0me[m,x]]] /. x -> R;

(* We don't have a known result to compare to, so instead can
integrate potential/density basis functions against each other to
check for orthogonality. E.g.,

   NIntegrate[R*sigmanme[n,m,R]*psinme[n,m,R]/Nnme[n,m], {R,0,∞}]

should give -1.0 for all (n,m) *)


(* These commands produce CSV files used in the plots in Fig. 2: *)

psitable = Table[FullSimplify[-psinme[n,0,r]/Sqrt[Nnme[n,0]], Assumptions -> {Element[r,Reals],r>0}], {n,0,3}];

psitable2 = Table[FullSimplify[-psinme[n,1,r]/Sqrt[Nnme[n,1]], Assumptions -> {Element[r,Reals],r>0}], {n,0,3}];

d = Table[{r, psitable[[1]],psitable[[2]],psitable[[3]],psitable[[4]]} /. r -> 10^t, {t,-3,3,0.05}]; Export["expdisk_potentials.dat", d, "CSV"];

d2 = Table[{r, psitable2[[1]],psitable2[[2]],psitable2[[3]],psitable2[[4]]} /. r -> 10^t, {t,-3,3,0.05}]; Export["expdisk_potentials_m1.dat", d2, "CSV"];


