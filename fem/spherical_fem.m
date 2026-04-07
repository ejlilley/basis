(* Validation for integrals of the form I_{knl}(a,b) = \int_a^b r^k \Phi_{nl}(r),
   for Zhao-type, Slater-type & Gaussian-type Coulomb resolutions *)

Import["../coulomb/common.m"];

Import["../basis/spherical_expressions.m"];

ΦnlZ[n_,l_,r_] := With[{α = 1/2, ν = 1}, phinlA[3,α,ν,n,l,r]];

κnlZ[n_,l_] := With[{α = 1/2, ν = 1}, With[{μ = α(1+2*l)}, α^n*I^n*(n!) * Gamma[μ]*Gamma[μ+ν]/(Gamma[n+μ]*Gamma[μ+ν+n])*kn[n,μ/2,μ/2+ν]]];

βnlZ[n_,l_] := If[n==0, 0, With[{α = 1/2, ν = 1}, With[{μ = α(1+2*l)}, α^(-2) * βn[n,μ/2,μ/2+ν]]]];

PnlZ[n_,l_] := With[{α = 1/2, ν = 1}, With[{μ = α(1+2*l)}, I^n*(n!) * Gamma[μ]*Gamma[μ+ν]/(Gamma[n+μ]*Gamma[μ+ν+n])*pn[n,μ/2,μ/2+ν,α*s]]];

(* FullSimplify[phinl[α,1,0,l,r], Assumptions -> {Element[l,Integers],l>=0,Element[k,Integers],k>=0,Element[r,Reals],r>0,Element[α,Reals],α>=1/2,Element[ν,Reals],ν>=0}]; *)

Ik0lZ[k_,l_,r_] := With[{α = 1/2}, (r^(1 + k + l)*Hypergeometric2F1[(1 + k + l)*α, (1 + 2*l)*α, 1 + (1 + k + l)*α, -r^α^(-1)])/(1 + k + l)];

IknlZ[k_,-1,l_,a_,b_] := 0; IknlZ[k_,0,l_,a_,b_] := Ik0lZ[k,l,b] - Ik0lZ[k,l,a];

IknlZ[k_,n_,l_,a_,b_] := I*(b^(k+1)*ΦnlZ[n-1,l,b] - a^(k+1)*ΦnlZ[n-1,l,a])/κnlZ[n-1,l] - I*(k+1/2)*IknlZ[k,n-1,l,a,b] - βnlZ[n-1,l]*IknlZ[k,n-2,l,a,b];

IknlZtest[k_,n_,l_,a_,b_] := 1/κnlZ[n,l]*Integrate[r^k*ΦnlZ[n,l,r], {r,a,b}, Assumptions -> {Element[a,Reals],Element[b,Reals],a>0,b>0,b>a}];

TestIknlZ[k_,n_,l_,a_,b_] := FullSimplify[FunctionExpand[IknlZtest[k,n,l,a,b]/IknlZ[k,n,l,a,b]], Assumptions -> {Element[a,Reals],Element[b,Reals],a>0,b>0,b>a}];

(* TestIknlZ[6,5,3,a,b]

Out[241]= 1 *)


Import["../coulomb/slater.m"];

(* FullSimplify[FunctionExpand[Integrate[r^k*ΦnlS[0,l,r], r]], Assumptions -> {Element[l,Integers],l>=0,Element[k,Integers],k>=0,Element[r,Reals],r>0}]

FullSimplify[FunctionExpand[Integrate[r^l*ΦnlS[0,l,r], r]], Assumptions -> {Element[l,Integers],l>=0,Element[k,Integers],k>=0,Element[r,Reals],r>0}] // InputForm *)

Ik0lS[k_,l_,r_] := If[k==l,(-4*Pi*r^(1 + 2*l)*(1 + r - E^r*r^2*ExpIntegralE[-2*(1 + l), r] + (E^r*(1 + 2*l)*r^2*HypergeometricPFQ[{3 + 2*l, 3 + 2*l}, {4 + 2*l, 4 + 2*l}, -r])/(3 + 2*l)^2))/(E^r*(1 + 2*l)^2),(-4*Pi*(((k - l)*r^(1 + k + l)*(1 + r))/E^r + (1 + 2*l)*Gamma[3 + k + l, r] + (1 + k + l)*r^(k - l)*(Gamma[3 + 2*l] - Gamma[3 + 2*l, r])))/((k - l)*(1 + k + l)*(1 + 2*l))];

IknlS[k_,-1,l_,a_,b_] := 0; IknlS[k_,0,l_,a_,b_] := Ik0lS[k,l,b] - Ik0lS[k,l,a];

IknlS[k_,n_,l_,a_,b_] := I*(b^(k+1)*ΦnlS[n-1,l,b] - a^(k+1)*ΦnlS[n-1,l,a])/κnlS[n-1,l] - I*(k+1/2)*IknlS[k,n-1,l,a,b] - βnlS[n-1,l]*IknlS[k,n-2,l,a,b];

IknlStest[k_,n_,l_,a_,b_] := 1/κnlS[n,l]*Integrate[r^k*ΦnlS[n,l,r], {r,a,b}, Assumptions -> {Element[a,Reals],Element[b,Reals],a>0,b>0,b>a}];

TestIknlS[k_,n_,l_,a_,b_] := FullSimplify[FunctionExpand[IknlStest[k,n,l,a,b]/IknlS[k,n,l,a,b]], Assumptions -> {Element[a,Reals],Element[b,Reals],a>0,b>0,b>a}];

(* TestIknlS[7,4,5,a,b]

Out[118]= 1 *)


Import["../coulomb/gaussian.m"];

Ik0lG[k_,l_,r_] := If[k==l,(-2*Pi*r^(1 + 2*l)*HypergeometricPFQ[{1/2 + l, 1/2 + l}, {3/2 + l, 3/2 + l}, -r^2])/(1 + 2*l)^2,-((Pi*(r^k*(Gamma[1/2 + l] - Gamma[1/2 + l, r^2]) + r^l*Gamma[(1 + k + l)/2, r^2]))/((k - l)*r^l))];

IknlG[k_,-1,l_,a_,b_] := 0; IknlG[k_,0,l_,a_,b_] := Ik0lG[k,l,b] - Ik0lG[k,l,a];

IknlG[k_,n_,l_,a_,b_] := I*(b^(k+1)*ΦnlG[n-1,l,b] - a^(k+1)*ΦnlG[n-1,l,a])/κnlG[n-1,l] - I*(k+1/2)*IknlG[k,n-1,l,a,b] - βnlG[n-1,l]*IknlG[k,n-2,l,a,b];

IknlGtest[k_,n_,l_,a_,b_] := 1/κnlG[n,l]*Integrate[r^k*ΦnlG[n,l,r], {r,a,b}, Assumptions -> {Element[a,Reals],Element[b,Reals],a>0,b>0,b>a}];

TestIknlG[k_,n_,l_,a_,b_] := FullSimplify[FunctionExpand[IknlGtest[k,n,l,a,b]/IknlG[k,n,l,a,b]], Assumptions -> {Element[a,Reals],Element[b,Reals],a>0,b>0,b>a}];

(* TestIknlG[6,5,4,a,b]

Out[138]= 1 *)


