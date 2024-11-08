Import["common.m"];

ρ0l[l_,r_] := r^l * Exp[-r^2];

Φ0l[l_,r_] := (-Pi)*Gamma[l+1/2,0,r^2]/r^(l+1);

ωl[s_] := 1/8*Gamma[1/4 + l/2 + I*s/2]*Gamma[1/4 + l/2 - I*s/2];

Pnl[n_,l_,s_] := Pmp[n, s/2, l/2 + 1/4, π/2];

Testρnl[n_,l_,r_] := PolyOp[(Pnl[n,l,#])&, ζ[x], ρ0l[l,x]] /. x -> r;

TestΦnl[n_,l_,r_] := PolyOp[(Pnl[n,l,#])&, ζ2[x], Φ0l[l,x]] /. x -> r;

Nnl[n_,l_] := 1/4*hmp[n,l/2+1/4,Pi/2];

Fnl[n_,l_] := ((-I)^n*Pochhammer[l+1/2,n])/((n!));

Gnl[n_,l_] := 4*Pi*n!*(-1)^n/((2l+1)*Pochhammer[l+3/2,n]);

ρnl[n_,l_,r_] := I^n * r^l * Exp[-r^2] * (LaguerreL[n,l+1/2,2*r^2] + Boole[n>0]*LaguerreL[n-1,l+1/2,2*r^2]);

Φnl[n_,l_,r_] := Fnl[n,l]*(Φ0l[l,r] + Sum[Gnl[j,l]*r^l*Exp[-r^2]*LaguerreL[j,l+1/2,2*r^2], {j,0,n-1}]);

CheckPoisson[n_,l_,r_] := FullSimplify[(1/(4*Pi)*(r^(-2)*D[r^2*D[Φnl[n,l,r],r],r] - l*(l+1)/r^2*Φnl[n,l,r]))/(ρnl[n,l,r]), Assumptions -> {Element[l,Integers],l≥0,Element[r,Reals],r>=0}];

CheckOrthog[n_,nn_,l_] := NIntegrate[r^2*ρnl[n,l,r]*Conjugate[Φnl[nn,l,r]], {r,0,∞}]/(-Nnl[n,l]);
