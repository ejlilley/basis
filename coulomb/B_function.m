Import["common.m"];

(* B-funs *)

Bnl[0,l_,r_] := 2^(-l)/(l!)*r^(l-1)*Exp[-r];

Bnl[n_,l_,r_] := 1/(2^(n+l)*(n+l)!)*r^l*Exp[-r]*θn[n-1,r];

(* Potential corresponding to B-funs *)

βnl[0,l_,r_] := Gamma[(1+2*l),0,r]/r^(l+1)*(-4*Pi)/(2^l*(l!));

βnl[1,l_,r_] := (-4*Pi/((1+2l)))/(2^(l+1)*(l+1)!) * (Gamma[2,r]*r^(l) + Gamma[3+2*l,0,r]/r^(l+1));

βnlgeneral[n_,l_,r_] := -Exp[-I*Pi*n]*Sec[Pi*n]*4*Pi*Sqrt[2/Pi]*(2^(-5/2 - l - 2*n)*Pi*r^(l)* (2^(1 + 2*n)*Gamma[1/2 + l]*HypergeometricPFQRegularized[{1/2 + l}, {3/2 + l, 1/2 - n}, r^2/4] - r^(1 + 2*n)*Gamma[1 + l + n]* HypergeometricPFQRegularized[{1 + l + n}, {3/2 + n, 2 + l + n}, r^2/4])* (-1)^n)/Gamma[1 + l + n];

(* depending on whether n is integer, half-integer, or neither, select the correct definition of βnl *)

βnl[n_,l_,r_] := If[IntegerQ[n],βnl[n-1,l,r] + (2^(2 - l - n)*Pi)/((l+n)!)*r^l*Exp[-r]*θn[n-1,r],If[IntegerQ[n+1/2],(-1)^(n+1/2)/Pi*D[βnlgeneral[m,l,r]*Cos[Pi*m],m] /. m -> n,βnlgeneral[n,l,r]]];

(* Mellin-space polynomial & weight fun. corresponding to B-fun *)

ωl[q_,l_,s_] := (Gamma[(1 + 2*l - (2*I)*s)/4]*Gamma[(1 + 2*l + (2*I)*s)/4]*Gamma[(3 + 2*l - (2*I)*s)/4 + q]*Gamma[(3 + 2*l + (2*I)*s)/4 + q])/(2*Pi*Gamma[1 + l + q]^2);

Pnl[q_,n_,l_,s_] := pn[n,l/2+3/4+q,l/2+1/4,s/2];

(* Φnl and ρnl as evaluated by explicitly applying Pnl to the zeroth-orders *)

Testρnl[q_,n_,l_,r_] := PolyOp[(Pnl[q,n,l,#])&, ζ[x], Bnl[q,l,x]] /. x -> r;

TestΦnl[q_,n_,l_,r_] := PolyOp[(Pnl[q,n,l,#])&, ζ2[x], βnl[q,l,x]] /. x -> r;

(* Φnl and ρnl as evaluated by taking a linear combination of Bnl or βnl functions *)

Testρnl2[q_,n_,l_,r_] := I^n*Pochhammer[l+3/2+2q,n]*Pochhammer[l+q+1,n]/(n!) * Sum[Pochhammer[-n,j]*Pochhammer[n+2l+2q+1,j]/((j!)*Pochhammer[l+2q+3/2,j]) * Bnl[j+q,l,r], {j,0,n}];

TestΦnl2[q_,n_,l_,r_] := I^n*Pochhammer[l+3/2+2q,n]*Pochhammer[l+q+1,n]/(n!) * Sum[Pochhammer[-n,j]*Pochhammer[n+2l+2q+1,j]/((j!)*Pochhammer[l+2q+3/2,j]) * βnl[j+q,l,r], {j,0,n}];

(* Limiting values of Φnl/r^l and ρnl/r^l as r -> 0 *)

ρnlLim[q_,n_,l_] :=  (I^n*2^(-2 - l)*(1 + 2*l^2 + q + 2*n*(1 + n + 2*q) + l*(3 + 4*n + 2*q))*Gamma[1/2 + l + n]*Gamma[-1/2 + q]*Gamma[1 + l + n + q])/(Sqrt[Pi]*n!*Gamma[3/2 + l]*Gamma[1 + l + q]*Gamma[2 + l + q]);

ΦnlLim[q_,n_,l_] := -((2^(1 - l)*I^n*Sqrt[Pi]*Gamma[1/2 + l + n]*Gamma[1/2 + q]*Gamma[1 + l + n + q])/(n!*Gamma[3/2 + l]*Gamma[1 + l + q]^2));

(* coeffs and other auxiliary quantities required for the explicit expressions for Φnl and ρnl *)

Jnmlq[q_,n_,m_,l_] := ((-1)^m*I^n*2^(-1 - l)*(l + 2*q)!*Gamma[q]*Gamma[1 + l + n + q]*Gamma[3 + 2*l + m + 2*q]*Gamma[3/2 + l + n + 2*q]*HypergeometricPFQRegularized[{-n, q, 1 + l + 2*q, 1 + 2*l + n + 2*q},{1 + l + q, -m + q, 3 + 2*l + m + 3*q}, 1])/(Sqrt[Pi]*m!*(l + q)!*Gamma[1 + n]);

Nnlq[q_,n_,l_] := Gamma[n+2l+2q+3]/(2^(2l+2q+3)*(n!));

NNnlq[q_,n_,l_] := 4/(2*Pi*Gamma[1 + l + q]^2) * hn[n,l/2+3/4+q,l/2+1/4];

Hnmlq[n_,j_,l_,q_] := (-1)^n*2^(q-1)*(l+q)!*Gamma[5/2+l+2q+n]*Gamma[3+j+2l+2q]/(Sqrt[Pi]*j!*(1+l+q)*(1+2l+2q)*Pochhammer[l+3/2,n]*Gamma[2+2l+n+2q])*Sum[(-1)^k*Binomial[n+1,k+1]*(k+l+2q)!*Gamma[3+k+2l+n+2q]/((3+2k+2l+4q)*Gamma[2+k+l+q]*Gamma[2+k+2l+3q])*HypergeometricPFQ[{1, k - n, 3 + k + 2*l + n + 2*q}, {2 + k, 5/2 + k + l + 2*q}, 1]*HypergeometricPFQ[{-j, 1 + 2*l + 2*q, 2 + 2*k + 2*l + 4*q},{3 + 2*l + 2*q, 2 + k + 2*l + 3*q}, 1], {k,0,n}];

(* the polynomial part of the Φnl expression *)

Γnlq[q_,0,l_,r_] := θn[q,r];

Γnlq[q_,n_,l_,r_] := Boole[EvenQ[n]]*Γnlq[q,0,l,r] + r^2*Sum[Boole[EvenQ[n-j]]*Sum[If[Or[m<0,j<0],0,(Hnmlq[j,m,l,q] - Hnmlq[j-2,m,l,q])/Nnlq[q,m,l]]*LaguerreL[m,2l+2q+2,2r], {m,j-q,j+q-2}], {j,1,n}];

(* compact forms of Φnl and ρnl *)

ρnl[q_,n_,l_,r_] := r^l*Exp[-r]*Sum[Jnmlq[q,n,j,l]/Nnlq[q,j,l]*LaguerreL[j,2l+2q+2,2*r], {j,Max[n-q-1,0],n+q-1}];

Φnl[q_,n_,l_,r_] := AAnlq[n,l,q]*(βnl[q,l,r] + Blq[q,l]*r^l*Exp[-r]*Γnlq[q,n-1,l,r]);

(* test functions (should evaluate to 1 for all reasonable combinations of parameters) *)

CheckPoisson[q_,n_,l_,r_] := FullSimplify[(r^(-2)*D[r^2*D[Φnl[q,n,l,r],r],r] - l*(l+1)/r^2*Φnl[q,n,l,r])/(4*Pi*ρnl[q,n,l,r])];

CheckOrthog[q_,n_,m_,l_] := -NIntegrate[FullSimplify[r^2*ρnl[q,n,l,r]*Conjugate[Φnl[q,m,l,r]]], {r,0,∞}]/N[NNnlq[q,n,l]];
