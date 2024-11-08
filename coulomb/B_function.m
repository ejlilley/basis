Import["common.m"];

(* B-funs *)

Bnl[0,l_,r_] := 2^(-l)/(l!)*r^(l-1)*Exp[-r];

Bnl[n_,l_,r_] := 1/(2^(n+l)*(n+l)!)*r^l*Exp[-r]*θn[n-1,r];

(* Potential corresponding to B-funs *)

βnl[0,l_,r_] := Gamma[(1+2*l),0,r]/r^(l+1)*(-4*Pi)/(2^l*(l!));

βnl[1,l_,r_] := (-4*Pi/((1+2l)))/(2^(l+1)*(l+1)!) * (Gamma[2,r]*r^(l) + Gamma[3+2*l,0,r]/r^(l+1));

βnl[n_,l_,r_] := -4*Pi*Sqrt[2/Pi]*(2^(-5/2 - l - 2*n)*Pi*r^(l)* (2^(1 + 2*n)*Gamma[1/2 + l]*HypergeometricPFQRegularized[{1/2 + l}, {3/2 + l, 1/2 - n}, r^2/4] - r^(1 + 2*n)*Gamma[1 + l + n]* HypergeometricPFQRegularized[{1 + l + n}, {3/2 + n, 2 + l + n}, r^2/4])* (-1)^n)/Gamma[1 + l + n];

(* Mellin-space polynomial & weight fun. corresponding to B-fun *)

ωl[q_,l_,s_] := (Gamma[(1 + 2*l - (2*I)*s)/4]*Gamma[(1 + 2*l + (2*I)*s)/4]*Gamma[(3 + 2*l - (2*I)*s)/4 + q]*Gamma[(3 + 2*l + (2*I)*s)/4 + q])/(2*Pi*Gamma[1 + l + q]^2);

Pnl[q_,n_,l_,s_] := pn[n,l/2+3/4+q,l/2+1/4,s/2];

Testρnl[q_,n_,l_,r_] := PolyOp[(Pnl[q,n,l,#])&, ζ[x], Bnl[q,l,x]] /. x -> r;

TestΦnl[q_,n_,l_,r_] := PolyOp[(Pnl[q,n,l,#])&, ζ2[x], βnl[q,l,x]] /. x -> r;

Jnmlq[q_,n_,m_,l_] := ((-1)^m*I^n*2^(-1 - l)*(l + 2*q)!*Gamma[q]*Gamma[1 + l + n + q]*Gamma[3 + 2*l + m + 2*q]*Gamma[3/2 + l + n + 2*q]*HypergeometricPFQRegularized[{-n, q, 1 + l + 2*q, 1 + 2*l + n + 2*q},{1 + l + q, -m + q, 3 + 2*l + m + 3*q}, 1])/(Sqrt[Pi]*m!*(l + q)!*Gamma[1 + n]);

Nnlq[q_,n_,l_] := Gamma[n+2l+2q+3]/(2^(2l+2q+3)*(n!));

NNnlq[q_,n_,l_] := 4/(2*Pi*Gamma[1 + l + q]^2) * hn[n,l/2+3/4+q,l/2+1/4];

AAnlq[n_,l_,q_] := (-I)^n*Pochhammer[l+1/2,n]*Pochhammer[1+l+q,n]/(n!);

Anl[n_,l_] := (-I)^n*Pochhammer[l+1/2,n]*Pochhammer[l+1,n]/(n!);

Blq[q_,l_] := 8*Pi/(2^(l+q)*((l+q)!)*(2l+1));

Cjl[j_,l_] := (-1)^j*(j!)/Pochhammer[2l+2,j];

Cnmlqu[n_,m_,l_,q_,u_] := (Sqrt[Pi]*Gamma[1 + l + m + q]*Gamma[3/2 + l + m + u]*Pochhammer[1/2 + l + m, -m + n]*Pochhammer[1 + l + m + q, -m + n]*Pochhammer[1 + 2*l + n + 2*q, m])/((-m + n)!*Gamma[(1 + m - n)/2]*Gamma[l + (2 + m + n)/2 + q]*Gamma[l + (3 + m + n)/2 + u]*Pochhammer[1 + 2*l + m + 2*u, m])*(Pochhammer[q-u,(n-m)/2]);

AAnjlq[n_,j_,l_,q_] := Sum[Cnmlqu[n,m,l,q,0]*Anl[m,l], {m,j,n}];

AAnjlq[n_,0,l_,q_] := AAnlq[n,l,q];

AAnjlq[n_,n_,l_,q_] := (-I)^n*Pochhammer[l+1/2,n]*Pochhammer[l+1,n]*Pochhammer[n+2l+2q+1,n]/(n!*Pochhammer[n+2l+1,n]);

Hnmlq[n_,m_,l_,q_] := If[Or[n<0,m<0],0,(2l+2q)!*(2q)!/(m!*q!*2^(2l+3q+1))*Sum[Pochhammer[-q,k]*Pochhammer[2l+2q+1,k]/(k!*Pochhammer[-2q,k])*Sum[AAnjlq[n+1,j+1,l,q]*Cjl[j,l]/AAnlq[n+1,l,q]*Pochhammer[2l+2,j]/(j!)*Sum[Pochhammer[-j,p]*Pochhammer[2l+2q+k+1,p]*Pochhammer[2-p-k,m]*Hypergeometric2F1[-p,k-q,k-2q,2]/((p!)*Pochhammer[2l+2,p]), {p,0,j}], {j,0,n}], {k,0,q}]];

Γnlq[q_,0,l_,r_] := θn[q,r];

Γnlq[q_,n_,l_,r_] := Boole[EvenQ[n]]*Γnlq[q,0,l,r] + r^2*Sum[Boole[EvenQ[n-j]]*Sum[If[Or[m<0,j<0],0,(Hnmlq[j,m,l,q] - Hnmlq[j-2,m,l,q])/Nnlq[q,m,l]]*LaguerreL[m,2l+2q+2,2r], {m,j-q,j+q-2}], {j,1,n}];

ρnl[q_,n_,l_,r_] := r^l*Exp[-r]*Sum[Jnmlq[q,n,j,l]/Nnlq[q,j,l]*LaguerreL[j,2l+2q+2,2*r], {j,Max[n-q-1,0],n+q-1}];

Φnl[q_,n_,l_,r_] := AAnlq[n,l,q]*(βnl[q,l,r] + Blq[q,l]*r^l*Exp[-r]*Γnlq[q,n-1,l,r]);

CheckPoisson[q_,n_,l_,r_] := FullSimplify[(r^(-2)*D[r^2*D[Φnl[q,n,l,r],r],r] - l*(l+1)/r^2*Φnl[q,n,l,r])/(4*Pi*ρnl[q,n,l,r])];

CheckOrthog[q_,n_,m_,l_] := -NIntegrate[FullSimplify[r^2*ρnl[q,n,l,r]*Conjugate[Φnl[q,m,l,r]]], {r,0,∞}]/N[NNnlq[q,n,l]];

