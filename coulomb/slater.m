Import["common.m"];

ρ0l[l_,r_] := r^l*Exp[-r];

Φ0l[l_,r_] := (-4*Pi/((1+2l))) * (Gamma[2,r]*r^(l) + Gamma[3+2*l,0,r]/r^(l+1));

ωl[l_,s_] := 2^(1+2l)/Pi*Gamma[1/4+l/2+I*s/2]*Gamma[1/4+l/2-I*s/2]*Gamma[7/4+l/2+I*s/2]*Gamma[7/4+l/2-I*s/2];

Pnl[n_,l_,s_] := pn[n,l/2+1/4,l/2+7/4,s/2];

Testρnl[n_,l_,r_] := PolyOp[(Pnl[n,l,#])&, ζ[x], ρ0l[l,x]] /. x -> r;

TestΦnl[n_,l_,r_] := PolyOp[(Pnl[n,l,#])&, ζ2[x], Φ0l[l,x]] /. x -> r;

Nnl[n_,l_] := 2^(3+2l)/Pi*hn[n,l/2+1/4,l/2+7/4];

AAnl[n_,l_] := I^n * Pochhammer[l+2,n] * Pochhammer[l+1/2,n]/((2l+1)*Pochhammer[2l+3,n]);

Snl[n_,l_] := (-I)^n*Pochhammer[l+1/2,n]*Pochhammer[l+2,n]/(n!);

Vml[m_,l_] := (2*Gamma[2l+3])*(2m+2l+5)*(m!)/Gamma[m+2l+5];

ρnl[n_,l_,r_] := AAnl[n,l]*r^l*Exp[-r]*((1+2l+2n)*LagL[n,2l+4,2r] - (5 + 2l + 2n)*LagL[n-2,2l+4,2r]);

Φnl[n_,l_,r_] := Snl[n,l]*(Φ0l[l,r] - 8*(-1)^n*Pi/(1+2l)*r^l*Exp[-r]*(Boole[OddQ[n]]*(1+r) - r^2*Sum[Boole[EvenQ[n-m]]*Vml[m,l]*LagL[m,2l+4,2*r],{m,0,n-2}]));

CheckPoisson[n_,l_,r_] := FullSimplify[(r^(-2)*D[r^2*D[Φnl[n,l,r],r],r] - l*(l+1)/r^2*Φnl[n,l,r])/(4*Pi*ρnl[n,l,r])];

CheckOrthog[n_,m_,l_] := -NIntegrate[r^2*ρnl[n,l,r]*Conjugate[Φnl[m,l,r]]/Nnl[n,l], {r,0,∞}];

