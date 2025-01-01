Import["common.m"];

(* Alternative (worse) method of computing the coefficients in Φnl *)

AAnlq[n_,l_,q_] := (-I)^n*Pochhammer[l+1/2,n]*Pochhammer[1+l+q,n]/(n!);

Anl[n_,l_] := (-I)^n*Pochhammer[l+1/2,n]*Pochhammer[l+1,n]/(n!);

Blq[q_,l_] := 8*Pi/(2^(l+q)*((l+q)!)*(2l+1));

Cjl[j_,l_] := (-1)^j*(j!)/Pochhammer[2l+2,j];

Cnmlqu[n_,m_,l_,q_,u_] := (Sqrt[Pi]*Gamma[1 + l + m + q]*Gamma[3/2 + l + m + u]*Pochhammer[1/2 + l + m, -m + n]*Pochhammer[1 + l + m + q, -m + n]*Pochhammer[1 + 2*l + n + 2*q, m])/((-m + n)!*Gamma[(1 + m - n)/2]*Gamma[l + (2 + m + n)/2 + q]*Gamma[l + (3 + m + n)/2 + u]*Pochhammer[1 + 2*l + m + 2*u, m])*(Pochhammer[q-u,(n-m)/2]);

AAnjlq[n_,j_,l_,q_] := Sum[Cnmlqu[n,m,l,q,0]*Anl[m,l], {m,j,n}];

AAnjlq[n_,0,l_,q_] := AAnlq[n,l,q];

AAnjlq[n_,n_,l_,q_] := (-I)^n*Pochhammer[l+1/2,n]*Pochhammer[l+1,n]*Pochhammer[n+2l+2q+1,n]/(n!*Pochhammer[n+2l+1,n]);

Hnmlq[n_,m_,l_,q_] := If[Or[n<0,m<0],0,(2l+2q)!*(2q)!/(m!*q!*2^(2l+3q+1))*Sum[Pochhammer[-q,k]*Pochhammer[2l+2q+1,k]/(k!*Pochhammer[-2q,k])*Sum[AAnjlq[n+1,j+1,l,q]*Cjl[j,l]/AAnlq[n+1,l,q]*Pochhammer[2l+2,j]/(j!)*Sum[Pochhammer[-j,p]*Pochhammer[2l+2q+k+1,p]*Pochhammer[2-p-k,m]*Hypergeometric2F1[-p,k-q,k-2q,2]/((p!)*Pochhammer[2l+2,p]), {p,0,j}], {j,0,n}], {k,0,q}]];

