Import["common.m"];

ρ0lS[l_,r_] := r^l*Exp[-r];

Φ0lS[l_,r_] := (-4*Pi/((1+2l))) * (Gamma[2,r]*r^(l) + Gamma[3+2*l,0,r]/r^(l+1));

ωlS[l_,s_] := 2^(1+2l)/Pi*Gamma[1/4+l/2+I*s/2]*Gamma[1/4+l/2-I*s/2]*Gamma[7/4+l/2+I*s/2]*Gamma[7/4+l/2-I*s/2];

PnlS[n_,l_,s_] := pn[n,l/2+1/4,l/2+7/4,s/2];

TestρnlS[n_,l_,r_] := PolyOp[(PnlS[n,l,#])&, ζ[x], ρ0lS[l,x]] /. x -> r;

TestΦnlS[n_,l_,r_] := PolyOp[(PnlS[n,l,#])&, ζ2[x], Φ0lS[l,x]] /. x -> r;

κnlS[n_,l_] := If[n<0,0,kn[n,l/2+1/4,l/2+7/4]/2^n];

βnlS[n_,l_] := 4*βn[n,l/2+1/4,l/2+7/4];

NnlS[n_,l_] := 2^(3+2l)/Pi*hn[n,l/2+1/4,l/2+7/4];

AAnl[n_,l_] := I^n * Pochhammer[l+2,n] * Pochhammer[l+1/2,n]/((2l+1)*Pochhammer[2l+3,n]);

Snl[n_,l_] := (-I)^n*Pochhammer[l+1/2,n]*Pochhammer[l+2,n]/(n!);

Vml[m_,l_] := (2*Gamma[2l+3])*(2m+2l+5)*(m!)/Gamma[m+2l+5];

ρnlS[n_,l_,r_] := AAnl[n,l]*r^l*Exp[-r]*((1+2l+2n)*LagL[n,2l+4,2r] - (5 + 2l + 2n)*LagL[n-2,2l+4,2r]);

ρnlmS[n_,l_,m_,r_,θ_,φ_] := If[n<0,0,ρnlS[n,l,r]*SphericalHarmonicY[l,m,θ,φ]];

ΦnlS[n_,l_,r_] := Snl[n,l]*(Φ0lS[l,r] - 8*(-1)^n*Pi/(1+2l)*r^l*Exp[-r]*(Boole[OddQ[n]]*(1+r) - r^2*Sum[Boole[EvenQ[n-m]]*Vml[m,l]*LagL[m,2l+4,2*r],{m,0,n-2}]));

CheckPoissonS[n_,l_,r_] := FullSimplify[(r^(-2)*D[r^2*D[ΦnlS[n,l,r],r],r] - l*(l+1)/r^2*ΦnlS[n,l,r])/(4*Pi*ρnlS[n,l,r])];

CheckOrthogS[n_,m_,l_] := -NIntegrate[r^2*ρnlS[n,l,r]*Conjugate[ΦnlS[m,l,r]]/NnlS[n,l], {r,0,∞}];

(* Linearization *)

CNjl[n_,j_,l_] := If[j>n,0,(-I)^j*Sqrt[Pi]*Gamma[2l+5]/2^(n+j+2l+4) * Pochhammer[-n,j]*Pochhammer[2l+3,j]*Pochhammer[2l+5,n]/(Pochhammer[l+3/2,j]*Pochhammer[l+2,j]) * 4/(2j+2l+5)*(D[z^((2j+2l+5)/4)*HypergeometricPFQRegularized[{1, (j - n)/2, (1 + j - n)/2}, {(5 + j)/2 + l, 3 + j/2 + l}, z],z] /. z -> 1)];

Fnjk[l_,λ_,ll_,n_,ν_,-1,k_] := 0; Fnjk[l_,λ_,ll_,n_,ν_,j_,-1] := 0;

Fnjk[l_,λ_,ll_,-1,ν_,j_,k_] := 0; Fnjk[l_,λ_,ll_,n_,-1,j_,k_] := 0;

Fnjk[l_,λ_,ll_,0,0,j_,0] := CNjl[l+λ-ll,j,ll];

Fnjk[l_,λ_,ll_,n_,ν_,j_,k_] := Fnjk[l,λ,ll,n,ν,j,k] = With[{nn = l+λ-ll}, Piecewise[{{0, nn < 0}, {0, k > n+ν}, {0, j > nn + n + ν}, {Fnjk[λ,l,ll,ν,n,j,n+ν-k], ν > n}}, -I*(k+λ)*Fnjk[l,λ,ll,n-1,ν,j,k] - I*(k-l-(n-1)-ν-7/2)*Fnjk[l,λ,ll,n-1,ν,j,k-1] + κnlS[j-1,ll]/κnlS[j,ll]*Fnjk[l,λ,ll,n-1,ν,j-1,k] + κnlS[j+1,ll]*βnlS[j+1,ll]/κnlS[j,ll]*Fnjk[l,λ,ll,n-1,ν,j+1,k] - βnlS[n-1,l]*(Fnjk[l,λ,ll,n-2,ν,j,k] + 2*Fnjk[l,λ,ll,n-2,ν,j,k-1] + Fnjk[l,λ,ll,n-2,ν,j,k-2])]];

TestΛNjLnνS[l_,λ_,j_,ll_,n_,ν_,a_,b_] := a^l*b^λ/(a+b)^(l+λ+n+ν) * Sum[Fnjk[l,λ,ll,n,ν,j,k]*a^(n+ν-k)*b^(k), {k,0,n+ν}];

ΛnlνλNLS[n_,l_,ν_,λ_,nn_,ll_,a_,b_] := (-1)^ν*κnlS[n,l]*κnlS[ν,λ]/(κnlS[0,l]*κnlS[0,λ])*a^l*b^λ/(a+b)^(l+λ+n+ν) * Sum[Fnjk[l,λ,ll,n,ν,nn,k]*a^(n+ν-k)*b^(k), {k,0,n+ν}];

CheckLinearizationS[n_,l_,m_,ν_,λ_,μ_] := FullSimplify[ρnlmS[n,l,m,a*r,θ,φ]*Conjugate[ρnlmS[ν,λ,μ,b*r,θ,φ]] - Sum[Quiet[JLM[l,m,λ,μ,ll]]*(Sum[ΛnlνλNLS[n,l,ν,λ,nn,ll,a,b]*ρnlmS[nn,ll,m-μ,(a+b)*r,θ,φ], {nn,0,l+λ-ll+n+ν}]), {ll,Abs[l-λ],l+λ}], Assumptions -> {Element[r,Reals],r>=0,Element[a,Reals],a>=0,Element[b,Reals],b>=0,Element[l,Integers],l>=0,Element[θ,Reals],Element[φ,Reals],0<=θ,π>=θ,0<=φ,2*π>=φ}];

(* Table[Table[CheckLinearization[n,l,m,ν,λ,μ], {m,-l,l}, {μ,-λ,λ}], {n,0,2}, {ν,0,2}, {l,0,2}, {λ,0,2}] *)

(* Multiplication coefficients *)

MnlS[n_,l_,a_] := -Pi*(-I)^n*(n+l+1)*Gamma[2n+2l+1]*a^(l-2)/(2^(2n-3)*(n!)*(a+1)^(2l+3))* ((a-1)/(a+1))^n * (1 + (2n+2l+3)*a + a^2);

