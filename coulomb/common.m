(* radial Laplacian operator, operating on f[r]*Y_{lm} *)

Lap[l_,r_] = (r^(1-3)*D[r^(3-1)*D[#,r],r] - (l+3-2)*l/r^2*#)&;

(* the operator referred to as \mathcal{D} in the text *)

ζ[r_] = (I*(r*D[#,{r,1}] + 5/2*#))&;

(* \mathcal{D}^*, the adjoint operator of \mathcal{D} *)

ζ2[r_] = (I*(r*D[#,{r,1}] + 1/2*#))&;

(* manipulating spherical and cartesian coordinates *)

(* ReduceRange[a_,x_] := If[x < 0, ReduceRange[a,x+a], If[x > a, ReduceRange[a,x-a], x]];

CartToSpher[a_] := Block[{x = a[[1]], y = a[[2]], z = a[[3]]}, If[a == {0.0,0.0,0.0}, {0.0,0.0,0.0}, {Sqrt[x^2 + y^2 + z^2],ReduceRange[Pi,ArcTan[z,Sqrt[x^2+y^2]]],ReduceRange[2*Pi,ArcTan[x,y]]}]]; *)

CartToSpher[a_] := Block[{xxxx = a[[1]], yyyy = a[[2]], zzzz = a[[3]]}, {Sqrt[xxxx^2 + yyyy^2 + zzzz^2],ArcTan[zzzz,Sqrt[xxxx^2+yyyy^2]],ArcTan[xxxx,yyyy]}];

SpherToCart[a_] := Block[{rrrr = a[[1]], θθθθ = a[[2]], φφφφ = a[[3]]}, {rrrr*Sin[θθθθ]*Cos[φφφφ], rrrr*Sin[θθθθ]*Sin[φφφφ], rrrr*Cos[θθθθ]}];

AddSpherical[a_,b_] := CartToSpher[SpherToCart[a] + SpherToCart[b]];

SubSpherical[a_,b_] := CartToSpher[SpherToCart[a] - SpherToCart[b]];


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

θn[n_,x_] := If[n<0,θn[-n-1,x]*x^(2n+1),(n!)/(-2)^n*LaguerreL[n,-2*n-1,2*x]];

(* linearization coefficients for spherical harmonics *)

GLM[l1_,m1_,l2_,m2_,l_] := Boole[EvenQ[l1+l2-l]]*Boole[l>=Abs[l1-l2]]*Boole[l<=l1+l2]*Boole[Abs[m1+m2]<=l] * Sqrt[(2l1+1)*(2l2+1)*(2l+1)/(4*Pi)] * ThreeJSymbol[{l1,m1},{l2,m2},{l,-m1-m2}]*ThreeJSymbol[{l1,0},{l2,0},{l,0}]*(-1)^(m1+m2);

(* linearization coefficients for (spherical harmonic)*(c.c. spherical harmonic) *)

JLM[l1_,m1_,l2_,m2_,l_] := JLM[l1,m1,l2,m2,l] = GLM[l1,m1,l2,-m2,l]*(-1)^(m2);

CheckYlmLinearization[l1_,m1_,l2_,m2_] := FullSimplify[SphericalHarmonicY[l1,m1,θ,φ]*Conjugate[SphericalHarmonicY[l2,m2,θ,φ]] - Sum[Quiet[JLM[l1,m1,l2,m2,l]]*SphericalHarmonicY[l,m1-m2,θ,φ], {l,Abs[l1-l2],l1+l2}], Assumptions -> {Element[θ,Reals],Element[φ,Reals],0<=θ,π>=θ,0<=φ,2*π>=φ}];

(* linearization coefficients for Jacobi polynomials *)

JacLinCjsn[j_,s_,n_,a_,b_] := Sum[ Pochhammer[-j,r] * Pochhammer[j+2s+a+b+1,r]/((r!) * Pochhammer[a+b+2+2n,r]) * HypergeometricPFQ[{-r, r+2s+1, s-n, s-n-b}, {s+1, a+s+1, 2s-2n-a-b}, 1], {r,0,j}];

JacLinHjn[j_,s_,n_,a_,b_] := 2^(a+b+1) * Gamma[a+1] * Gamma[b+1] / Gamma[a+b+2+2n] * Pochhammer[b+1,n] * Pochhammer[a+b+1+n-s,n-s] * Pochhammer[a+b+1+s+j,s] * ((s+j)!) * (n!) / ( Pochhammer[a+1,s] * Pochhammer[a+1, n-s] * (s!) * (j!) ) * JacLinCjsn[j,s,n,a,b];

JacLinGkmn[k_,m_,n_,a_,b_] := 2^(-a-b-1) * Gamma[a+b+1] * Pochhammer[a+1,k] * Pochhammer[a+b+1,k] * (2k + a + b + 1) / (Gamma[a+1] * Gamma[b+1] * (k!) * Pochhammer[b+1,k]) * JacLinHjn[k-n+m,n-m,n,a,b];

JacLinLnmk[n_, m_, k_, a_, b_] := Pochhammer[a+1,n]/n! * Pochhammer[a+1,m]/m! * k!/Pochhammer[a+1,k] * If[m > n, JacLinGkmn[k,n,m,a,b], JacLinGkmn[k,m,n,a,b]];

CheckJacLin[n_,m_,a_,b_] := FullSimplify[JacobiP[n,a,b,x]*JacobiP[m,a,b,x] - Sum[ JacLinLnmk[n,m,k,a,b]*JacobiP[k,a,b,x], {k,Abs[n-m],n+m} ]];
