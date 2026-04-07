(* radial Laplacian operator, operating on f[r]*Y_{lm} *)

Lapl[l_,r_] = (r^(-2)*D[r^(2)*D[#,r],r] - (l+1)*l/r^2*#)&;

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

PolyOp[p_,op_,f_] := Block[{ss,jj}, Sum[Coefficient[p[ss], ss, jj]*Nest[op,f,jj], {jj,0,Exponent[p[ss],ss]}]];

(* Meixner-Pollaczek polynomial *)

Pmp[n_, x_, λ_, φ_] := Pochhammer[2*λ,n]/n! * Exp[I*n*φ] * Hypergeometric2F1[-n, λ + I*x, 2*λ, 1 - Exp[-2*I*φ]];

hmp[n_, λ_, φ_] := 2*π*Gamma[n+2*λ]/((2*Sin[φ])^(2*λ) * n!); (* normalisation constant *)

kmp[n_, λ_, φ_] := (2*Sin[φ])^n/(n!); (* coeff of x^n *)

(* linearization coefficients for symmetric Meixner-Pollaczek polynomials (i.e. M-P polys with φ=π/2), from Araaya (2004) *)

Lpqp[p_, q_, ν_, λ_] := Boole[EvenQ[p+q+ν]]*(ν!)*Gamma[2*λ + (p+q+ν)/2]/( Gamma[2*λ + ν] * Gamma[1 + (ν+p-q)/2] * Gamma[1 + (p+q-ν)/2] * Gamma[1 + (ν+q-p)/2] );

CheckMPLin[n_, m_, a_] := FullSimplify[Pmp[n, x, a, Pi/2]*Pmp[m, x, a, Pi/2] - Sum[Lpqp[n,m,j,a]*Pmp[j, x, a, Pi/2], {j, Abs[n-m], n+m}]];

(* symmetric continuous Hahn polynomial *)

pn[n_,a_,b_,x_] := I^n * Pochhammer[2*a,n]*Pochhammer[a+b,n]/(n!) * HypergeometricPFQ[{-n,n+2a+2b-1,a+I*x},{2a,a+b},1];

hn[n_,a_,b_] := 2*Pi*Gamma[n+2*a]*Gamma[n+2*b]*(Gamma[n+a+b])^2/((2*n+2*a+2*b-1)*Gamma[n+2*a+2*b-1]*(n!)); (* normalisation const for continuous Hahn polys *)

kn[n_,a_,b_] := Pochhammer[n + 2*a + 2*b - 1,n]/(n!);

βn[n_,a_,b_] := n*(n+2a-1)*(n+2b-1)*(n+2a+2b-2)/(4*(2n+2a+2b-3)*(2n+2a+2b-1));

(* Wilson polynomial *)

Wn[n_,a_,b_,c_,d_,x_] := Pochhammer[a+b,n]*Pochhammer[a+c,n]*Pochhammer[a+d,n]*HypergeometricPFQ[{-n,n+a+b+c+d-1,a+I*x,a-I*x},{a+b,a+c,a+d},1];

hnW[n_,a_,b_,c_,d_] := (n!)*2*Pi*Gamma[n+a+b]*Gamma[n+a+c]*Gamma[n+b+c]*Gamma[n+b+d]*Gamma[n+c+d]*Gamma[n+a+d]/((2*n-1+a+b+c+d)*Gamma[n-1+a+b+c+d]);

wW[a_,b_,c_,d_,x_] := (Gamma[a+I*x]*Gamma[b+I*x]*Gamma[c+I*x]*Gamma[d+I*x]/Gamma[2*I*x])*(Gamma[a-I*x]*Gamma[b-I*x]*Gamma[c-I*x]*Gamma[d-I*x]/Gamma[-2*I*x]);

(* Laguerre polynomial with non-negative degree *)

LagL[n_,a_,x_] := If[n < 0, 0, LaguerreL[n,a,x]];

(* reverse Bessel polynomial *)

θn[n_,x_] := If[n<0,θn[-n-1,x]*x^(2n+1),(n!)/(-2)^n*LaguerreL[n,-2*n-1,2*x]];

(* linearization coefficients for reverse Bessel polynomials, from arxiv:1109.4660 *)

Lθnmk[n_,m_,k_,a_,b_] := 2^(n+m-k) * Pochhammer[1/2,n+m-k] * a^(2n+m-k) * b^(k-n) * (a+b)^(-n-m) * Sum[(-1)^j*Binomial[n+m+1,2n+2m-2k-j]*Binomial[-m+k+j,j]*(1+b/a)^j, {j,0,2*(n+m-k)}];

(* test with e.g.: FullSimplify[θn[3,a*x]*θn[4,b*x] - Sum[Lθnmk[3,4,k,a,b]*θn[k,(a+b)*x], {k,0,3+4}]] *)

(* linearization coefficients for spherical harmonics *)

GLM[l1_,m1_,l2_,m2_,l_] := Boole[EvenQ[l1+l2-l]]*Boole[l>=Abs[l1-l2]]*Boole[l<=l1+l2]*Boole[Abs[m1+m2]<=l] * Sqrt[(2l1+1)*(2l2+1)*(2l+1)/(4*Pi)] * ThreeJSymbol[{l1,m1},{l2,m2},{l,-m1-m2}]*ThreeJSymbol[{l1,0},{l2,0},{l,0}]*(-1)^(m1+m2);

(* linearization coefficients for (spherical harmonic)*(c.c. spherical harmonic) *)

JLM[l1_,m1_,l2_,m2_,l_] := JLM[l1,m1,l2,m2,l] = GLM[l1,m1,l2,-m2,l]*(-1)^(m2);

CheckYlmLinearization[l1_,m1_,l2_,m2_] := FullSimplify[SphericalHarmonicY[l1,m1,θ,φ]*Conjugate[SphericalHarmonicY[l2,m2,θ,φ]] - Sum[Quiet[JLM[l1,m1,l2,m2,l]]*SphericalHarmonicY[l,m1-m2,θ,φ], {l,Abs[l1-l2],l1+l2}], Assumptions -> {Element[θ,Reals],Element[φ,Reals],0<=θ,π>=θ,0<=φ,2*π>=φ}];

Jlm[l_,m_] := 1/((2l+1)/(4*Pi) * (l-m)!/(l+m)!); (* sph harm normalisation constant *)

(* linearization coefficients for Jacobi polynomials *)

JacLinCjsn[j_,s_,n_,a_,b_] := JacLinCjsn[j,s,n,a,b] = Sum[ Pochhammer[-j,r] * Pochhammer[j+2s+a+b+1,r]/((r!) * Pochhammer[a+b+2+2n,r]) * HypergeometricPFQ[{-r, r+2s+1, s-n, s-n-b}, {s+1, a+s+1, 2s-2n-a-b}, 1], {r,0,j}];

JacLinHjn[j_,s_,n_,a_,b_] := JacLinHjn[j,s,n,a,b] = 2^(a+b+1) * Gamma[a+1] * Gamma[b+1] / Gamma[a+b+2+2n] * Pochhammer[b+1,n] * Pochhammer[a+b+1+n-s,n-s] * Pochhammer[a+b+1+s+j,s] * ((s+j)!) * (n!) / ( Pochhammer[a+1,s] * Pochhammer[a+1, n-s] * (s!) * (j!) ) * JacLinCjsn[j,s,n,a,b];

JacLinGkmn[k_,m_,n_,a_,b_] := JacLinGkmn[k,m,n,a,b] = 2^(-a-b-1) * Gamma[a+b+1] * Pochhammer[a+1,k] * Pochhammer[a+b+1,k] * (2k + a + b + 1) / (Gamma[a+1] * Gamma[b+1] * (k!) * Pochhammer[b+1,k]) * JacLinHjn[k-n+m,n-m,n,a,b];

JacLinLnmk[n_, m_, k_, a_, b_] := Pochhammer[a+1,n]/n! * Pochhammer[a+1,m]/m! * k!/Pochhammer[a+1,k] * If[m > n, JacLinGkmn[k,n,m,a,b], JacLinGkmn[k,m,n,a,b]];

CheckJacLin[n_,m_,a_,b_] := FullSimplify[JacobiP[n,a,b,x]*JacobiP[m,a,b,x] - Sum[ JacLinLnmk[n,m,k,a,b]*JacobiP[k,a,b,x], {k,Abs[n-m],n+m} ]];


(* fast compiled Laguerre polynomial evaluation *)

genlaguerreC[n_,a_] := genlaguerreC[n,a] = Compile[{{x, _Real}}, Module[{kk, p, d, k, bi,i}, If[a<-1,1.0/0.0,If[n<0,0.0,If[n==0,1.0,If[n==1,-x+a+1,(d=-x/(a+1); p=d+1; Do[k=kk+1; d=-x/(k+a+1)*p + (k/(k+a+1))*d; p=d+p, {kk,0,n-2}]; bi = 1.0; Do[bi = bi*(n+a-i+1)/i, {i,1,n}]; bi*p)]]]]], RuntimeOptions -> {"EvaluateSymbolically" -> False}, RuntimeAttributes -> {Listable}, Parallelization -> True, CompilationTarget -> "C", CompilationOptions -> {"InlineCompiledFunctions" -> False, "InlineExternalDefinitions" -> True}];

(* same but output Laguerre of deg n and n-1 at same time *)

genlaguerreC2[n_,a_] := genlaguerreC2[n,a] = Compile[{{x, _Real}}, Module[{kk, p, p2, d, k, bi, bi2, i}, If[a<-1,{1.0/0.0,1.0/0.0},If[n<0,{0.0,0.0},If[n==0,{1.0,0.0},If[n==1,{-x+a+1,1.0},(d=-x/(a+1); p2=0.0; p=d+1; Do[k=kk+1; d=-x/(k+a+1)*p + (k/(k+a+1))*d; p2=p; p=d+p, {kk,0,n-2}]; bi=1.0; Do[bi = bi*(n+a-i+1)/i, {i,1,n}]; bi2=n/(n+a)*bi; {bi*p, bi2*p2})]]]]], RuntimeOptions -> {"EvaluateSymbolically" -> False}, RuntimeAttributes -> {Listable}, Parallelization -> True, CompilationTarget -> "C", CompilationOptions -> {"InlineCompiledFunctions" -> False, "InlineExternalDefinitions" -> True}];

(* same but output first n Laguerre polynomials *)

(* genlaguerreCN[n_,a_] := genlaguerreCN[n,a] = Compile[{{x, _Real}}, Module[{kk, p, d, k, bi, i, bi2={1.0}, p2={1.0}}, If[a<-1,{1.0/0.0},If[n<0,{},If[n==0,{1.0},If[n==1,{-x+a+1,1.0},(d=-x/(a+1); p=d+1; AppendTo[p2, p]; Do[k=kk+1; d=-x/(k+a+1)*p + (k/(k+a+1))*d; p=d+p; AppendTo[p2,p], {kk,0,n-2}]; bi=1.0; Do[bi = bi*(a+i)/i; AppendTo[bi2,bi], {i,1,n}]; Times[bi2 * p2])]]]]], RuntimeOptions -> {"EvaluateSymbolically" -> False}, RuntimeAttributes -> {Listable}, Parallelization -> True, CompilationTarget -> "C", CompilationOptions -> {"InlineCompiledFunctions" -> False, "InlineExternalDefinitions" -> True}]; *)

(* genlaguerreCN = Compile[{{n, _Integer},{a, _Real},{x, _Real}}, Module[{kk, p, d, k, bi, i, bi2={1.0}, p2={1.0}}, If[a<-1,{1.0/0.0},If[n<0,{},If[n==0,{1.0},If[n==1,{1.0,-x+a+1},(d=-x/(a+1); p=d+1; AppendTo[p2, p]; Do[k=kk+1; d=-x/(k+a+1)*p + (k/(k+a+1))*d; p=d+p; AppendTo[p2,p], {kk,0,n-2}]; bi=1.0; Do[bi = bi*(a+i)/i; AppendTo[bi2,bi], {i,1,n}]; Times[bi2 * p2])]]]]], RuntimeOptions -> {"EvaluateSymbolically" -> False}, RuntimeAttributes -> {Listable}, Parallelization -> True, CompilationTarget -> "C", CompilationOptions -> {"InlineCompiledFunctions" -> False, "InlineExternalDefinitions" -> True}]; *)

genlaguerreCN = Compile[{{n, _Integer},{a, _Real},{x, _Real}}, Module[{kk, p, d, k, bi, i, bi2={1.0}, p2={1.0}}, If[a<-1,{1.0/0.0},If[n<0,{},If[n==0,{1.0},If[n==1,{1.0,-x+a+1},(d=-x/(a+1); p=d+1; AppendTo[p2, p]; Do[k=kk+1; d=-x/(k+a+1)*p + (k/(k+a+1))*d; p=d+p; AppendTo[p2,p], {kk,0,n-2}]; bi=1.0; Do[bi = bi*(a+i)/i; AppendTo[bi2,bi], {i,1,n}]; Times[bi2 * p2])]]]]], RuntimeOptions -> {"EvaluateSymbolically" -> False}, RuntimeAttributes -> {Listable}, Parallelization -> True, CompilationOptions -> {"InlineCompiledFunctions" -> False, "InlineExternalDefinitions" -> True}];


(* Apply Im or Re to the second element in a list of pairs *)

Im2[s_] := Map[({#[[1]], Im[#[[2]]]})&, s];

Re2[s_] := Map[({#[[1]], Re[#[[2]]]})&, s];

scalingFun[y_] := Function[x, If[NumberQ[Log[Abs[x-1]]], Log[10,Abs[x-1]], y]];


