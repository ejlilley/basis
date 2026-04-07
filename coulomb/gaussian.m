Import["common.m"];

ρ0lG[l_,r_] := r^l * Exp[-r^2];

Φ0lG[l_,r_] := (-Pi)*Gamma[l+1/2,0,r^2]/r^(l+1);

ωlG[s_] := 1/8*Gamma[1/4 + l/2 + I*s/2]*Gamma[1/4 + l/2 - I*s/2];

PnlG[n_,l_,s_] := Pmp[n, s/2, l/2 + 1/4, π/2];

TestρnlG[n_,l_,r_] := PolyOp[(PnlG[n,l,#])&, ζ[x], ρ0lG[l,x]] /. x -> r;

TestΦnlG[n_,l_,r_] := PolyOp[(PnlG[n,l,#])&, ζ2[x], Φ0lG[l,x]] /. x -> r;

βnlG[n_,l_] := n(n + l − 1/2);

κnlG[n_,l_] := If[n<0,1,kmp[n,l/2+1/4,Pi/2]/2^n];

NnlG[n_,l_] := 1/4*hmp[n,l/2+1/4,Pi/2];

Fnl[n_,l_] := ((-I)^n*Pochhammer[l+1/2,n])/((n!));

Gnl[n_,l_] := 4*Pi*n!*(-1)^n/((2l+1)*Pochhammer[l+3/2,n]);

ρnlG[n_,l_,r_] := I^n * r^l * Exp[-r^2] * (LaguerreL[n,l+1/2,2*r^2] + Boole[n>0]*LaguerreL[n-1,l+1/2,2*r^2]);

ρnlmG[n_,l_,m_,r_,θ_,φ_] := If[n<0,0,ρnlG[n,l,r]*SphericalHarmonicY[l,m,θ,φ]];

ΦnlG[n_,l_,r_] := Fnl[n,l]*(Φ0lG[l,r] + Sum[Gnl[j,l]*r^l*Exp[-r^2]*LaguerreL[j,l+1/2,2*r^2], {j,0,n-1}]);

ΦnlmG[n_,l_,m_,r_,θ_,φ_] := If[n<0,0,ΦnlG[n,l,r]*SphericalHarmonicY[l,m,θ,φ]];

CheckPoissonG[n_,l_,r_] := FullSimplify[(1/(4*Pi)*(r^(-2)*D[r^2*D[ΦnlG[n,l,r],r],r] - l*(l+1)/r^2*ΦnlG[n,l,r]))/(ρnlG[n,l,r]), Assumptions -> {Element[l,Integers],l≥0,Element[r,Reals],r>=0}];

CheckOrthogG[n_,nn_,l_] := NIntegrate[r^2*ρnlG[n,l,r]*Conjugate[ΦnlG[nn,l,r]], {r,0,∞}]/(-NnlG[n,l]);


(* Linearization *)

Anjl[n_,j_,l_] := Pochhammer[-n,j]*Gamma[n+l+3/2]/(2^n*(I)^j*Gamma[j+l+3/2])*Hypergeometric2F1[1,j-n,j+l+3/2,-1];

Fnjk[l_,λ_,ll_,n_,ν_,-1,k_] := 0; Fnjk[l_,λ_,ll_,n_,ν_,j_,-1] := 0; Fnjk[l_,λ_,ll_,-1,ν_,j_,k_] := 0; Fnjk[l_,λ_,ll_,n_,-1,j_,k_] := 0;

Fnjk[l_,λ_,ll_,0,0,j_,0] := Anjl[(l+λ-ll)/2,j,ll];

Fnjk[l_,λ_,ll_,n_,ν_,j_,k_] := Fnjk[l,λ,ll,n,ν,j,k] = With[{nn = (l+λ-ll)/2}, Piecewise[{{0, Not[IntegerQ[nn]]}, {0, nn < 0}, {0, k > n+ν}, {0, j > nn + n + ν}, {Fnjk[λ,l,ll,ν,n,j,n+ν-k], ν > n}}, -I*(λ+2k)*Fnjk[l,λ,ll,n-1,ν,j,k] - I/2*(4k-2l-4n-4ν-5)*Fnjk[l,λ,ll,n-1,ν,j,k-1] - βnlG[n-1,l]*Fnjk[l,λ,ll,n-2,ν,j,k] - 2*βnlG[n-1,l]*Fnjk[l,λ,ll,n-2,ν,j,k-1] - βnlG[n-1,l]*Fnjk[l,λ,ll,n-2,ν,j,k-2] + κnlG[j-1,ll]/κnlG[j,ll]*Fnjk[l,λ,ll,n-1,ν,j-1,k] + κnlG[j+1,ll]*βnlG[j+1,ll]/κnlG[j,ll]*Fnjk[l,λ,ll,n-1,ν,j+1,k]]];

(* the linearization coefficients (note the factor of (-1)^ν appears because we are linearising ρ_{nlm} * Conjugate[ρ_{νλμ}] rather than just ρ_{nlm} * ρ_{νλμ} *)

ΛnlνλNLG[n_,l_,ν_,λ_,nn_,ll_,a_,b_] := (-1)^ν*κnlG[n,l]*κnlG[ν,λ]/(κnlG[0,l]*κnlG[0,λ]) * a^l*b^λ/(Sqrt[a^2+b^2])^(l+λ+2n+2ν) * Sum[Fnjk[l,λ,ll,n,ν,nn,k]*a^(2n+2ν-2k)*b^(2k), {k,0,n+ν}];

CheckLinearizationG[n_,l_,m_,ν_,λ_,μ_] := FullSimplify[ρnlmG[n,l,m,a*r,θ,φ]*Conjugate[ρnlmG[ν,λ,μ,b*r,θ,φ]] - Sum[Quiet[JLM[l,m,λ,μ,ll]]*If[OddQ[l+λ-ll],0,(Sum[ΛnlνλNLG[n,l,ν,λ,nn,ll,a,b]*ρnlmG[nn,ll,m-μ,Sqrt[a^2+b^2]*r,θ,φ], {nn,0,(l+λ-ll)/2+n+ν}])], {ll,Abs[l-λ],l+λ}], Assumptions -> {Element[r,Reals],r>=0,Element[a,Reals],a>=0,Element[b,Reals],b>=0,Element[l,Integers],l>=0,Element[θ,Reals],Element[φ,Reals],0<=θ,π>=θ,0<=φ,2*π>=φ}];


(* test it by running e.g.
   Table[Table[CheckLinearization[n,l,m,ν,λ,μ], {m,-l,l}, {μ,-λ,λ}], {n,0,3}, {ν,0,3}, {l,0,3}, {λ,0,3}]
   and go and make a coffee :-) *)


(* Addition *)

(* constant in Fourier transform of density basis function *)

Anlt[n_,l_] := Pi^(3/2)/(2^l * I^(n+l));

(* connection coefficient that doubles the argument of a Laguerre polynomial *)

djNL[j_,nn_,ll_] := Pochhammer[ll+1/2,j] * nn! / ( 2^j * j! * Pochhammer[ll+1/2,nn] ) * Binomial[j,nn];

(* linearization coefficient between 3 general Laguerre polynomials of parameters (l-1/2,λ-1/2,L-1/2) *)

Lnmk[n_,l_,ν_,λ_,k_,ll_] := (-1)^(n+ν-k)*(k!)/(ν!*n!) * Sum[ Binomial[n,j] * Binomial[ν,k-j] * Sum[Binomial[n-j,p] * Binomial[ν+j-k,p] * (p!) * Pochhammer[k+(ll-1/2)+1,p] * Pochhammer[ll-l+k-n+p+1,n-j-p] * Pochhammer[ll-λ+k-ν+p+1,ν+j-k-p], {p,0,Min[n-j,ν+j-k]}], {j,Max[0,k-ν], Min[k,n]}];

(* radial coefficient in the formula for addition theorem *)

CnlνλNL[n_,l_,ν_,λ_,nn_,ll_] := Anlt[n,l]*Conjugate[Anlt[ν,λ]]/(2^((ll+1)/2)*Anlt[nn,ll]) * Sum[djNL[j,nn,ll]*Lnmk[n,l,ν,λ,j,ll], {j,nn,n+ν}];

(* the addition theorem *)

Anlmνλμ[n_,l_,m_,ν_,λ_,μ_,r_,θ_,φ_] := -FullSimplify[Sum[Quiet[JLM[l,m,λ,μ,ll]]*If[OddQ[l+λ-ll],0,Sum[CnlνλNL[n,l,ν,λ,nn,ll]*(-1)^((l+λ-ll)/2)*(If[l+λ-ll == 0, ΦnlG[nn,ll,rd/Sqrt[2]], 2*Pi*Nest[Lapl[ll,rd],ρnlG[nn,ll,rd/Sqrt[2]],(l+λ-ll)/2 - 1]] /. rd -> r)*SphericalHarmonicY[ll,m-μ,θ,φ],{nn,0,n+ν}]],{ll,Abs[l-λ],l+λ}], Assumptions -> {Element[r,Reals],r>0,Element[θ,Reals],θ>=0,θ<=Pi,Element[φ,Reals],φ>=0,φ<=2*Pi}];

(* the same quantity as a numerical integral *)

TestAnlmνλμ[n_,l_,m_,ν_,λ_,μ_,r_,θ_,φ_] := -With[{xx = r*Sin[θ]*Cos[φ], yy = r*Sin[θ]*Sin[φ], zz = r*Cos[θ]}, Quiet[NIntegrate[ρnlm[n,l,m,Sqrt[(x+xx)^2 + (y+yy)^2 + (z+zz)^2],ArcTan[z+zz,Sqrt[(x+xx)^2+(y+yy)^2]],ArcTan[x+xx,y+yy]]*Conjugate[Φnlm[ν,λ,μ,Sqrt[x^2 + y^2 + z^2],ArcTan[z,Sqrt[x^2+y^2]],ArcTan[x,y]]], {x,-∞,∞}, {y,-∞,∞}, {z,-∞,∞}]]];

(* Test it!

   Clear[n,l,m,ν,λ,μ,r,θ,φ];
   
   n = 4; l = 2; m = -2; ν = 3; λ = 5; μ = -4; r = 2.193; θ = 0.318; φ = 1.849;
   
   Anlmνλμ[n,l,m,ν,λ,μ,r,θ,φ]
   
   TestAnlmνλμ[n,l,m,ν,λ,μ,r,θ,φ]

   *)

(* Multiplication coefficients *)

MnlG[n_,l_,a_] := -(I^n)*Gamma[n+l+1/2]*a^(l-2)/(2*(n!)*(a^2+1)^(l+1/2))*((a^2-1)/(a^2+1))^n;

MnνlG[n_,ν_,l_,a_] := Pi/2*I^(n+ν)*Gamma[n+ν+l+1/2]/((n!)*(ν!))*a^(l-2)/(a^2+1)^(l+1/2)*(((-1 + a^2)/(1 + a^2)))^(n+ν)*Hypergeometric2F1[-n,-ν,1/2-l-n-ν,((-1 + a^2)^2/(1 + a^2)^2)^(-1)];
