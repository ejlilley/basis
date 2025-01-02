(* Mathematica code to produce symbolic expressions for the spherical-coordinate
basis sets mentioned in the paper *)

ζ[r_] = (I*(r*D[#,{r,1}] + 5/2*#))&; (* The operator referred to as "D" in the paper (Eq. (2.4)) *)

ζ2[r_] = (I*(r*D[#,{r,1}] + 1/2*#))&; (* same operator but minus 2i *)

(* Symmetric continuous Hahn polynomials (see App. E.2) *)

pn[n_,a_,b_,x_] := I^n * Pochhammer[2*a,n]*Pochhammer[a+b,n]/(n!) * HypergeometricPFQ[{-n,n+2a+2b-1,a+I*x},{2a,a+b},1];

hn[n_,a_,b_] := 2*Pi*Gamma[n+2*a]*Gamma[n+2*b]*(Gamma[n+a+b])^2/((2*n+2*a+2*b-1)*Gamma[n+2*a+2*b-1]*(n!)); (* normalisation const for continuous Hahn polys *)

(***********************************************************)

(* Clutton-Brock (1973)'s basis set (see Sec. 4.1.1) *)

Phi0lCB73[l_,r_] := -r^l/(1+r^2)^(l+1/2);

rho0lCB73[l_,r_] := (1+2l)*(3+2l)/(4*Pi)*r^l/(1+r^2)^(l+5/2);

KnlCB73[n_,l_] := (1 + 2 l + 2 n) (3 + 2 l + 2 n); (* equal to 4n(n+2l+2) + (2l+1)(2l+3) *)

NnlCB73[n_,l_] := (n+2l+1)!/(2^(4l+6)*(n+l+1)*(n!)*((l)!)^2);

PhinlCB73[n_,l_,r_] := -r^l/(1 + r^2)^(l+1/2) * GegenbauerC[n,l+1,(r^2-1)/(r^2+1)];

rhonlCB73[n_,l_,r_] := -KnlCB73[n,l]/(4*Pi) * PhinlCB73[n,l,r]/(1+r^2)^2;

ωlCB73[l_, s_] := Re[((l+1/2)^2 + s^2)/(8*Pi^2) * MellinTransform[Phi0lCB73[l,r], r, 1/2 + I*s] * MellinTransform[Phi0lCB73[l,r], r, 1/2 - I*s]]; (* polynomial weight-function Eq (4.3) *)

PnlCB73[n_,l_,s_] := I^(n) * Sqrt[Pi] * Gamma[l+1/2] * (n+2l+1)!/(2^(2l) * (2n+2l+1)*(l)!*Gamma[n+l+1/2]^2) * pn[n,1/4 + l/2,5/4 + l/2,s/2]; (* the index-raising polynomial (see Eq. (4.4)); the complicated prefactor is just for compatibility with CB73's original expressions *)

PζrhoCB73[n_,l_,r_] := FullSimplify[Sum[Coefficient[PnlCB73[n,l,s], s, j]*Nest[ζ[x],rho0lCB73[l,x],j], {j,0,n}] /. x -> r]; (* the density basis function as expressed explicitly via the index-raising polynomial *)

PζPhiCB73[n_,l_,r_] := FullSimplify[Sum[Coefficient[PnlCB73[n,l,s], s, j]*Nest[ζ2[x],Phi0lCB73[l,x],j], {j,0,n}] /. x -> r];

(* Now can check the result by evaluating e.g.
   FullSimplify[PζPhiCB73[n,l,r]/PhinlCB73[n,l,r]]
   for various values of (n,l)
*)

(***********************************************************)

(* Lilley, Sanders & Evans (2018)'s 'A' basis set (Sec. (4.1.2) but
see the original paper for the basis functions themselves *)

phinlA[2,α_,ν_,0,0,r_] := 1/ν * (2ν-1) * (1 + r^(1/α))^(-ν) * Hypergeometric2F1[ν,1,ν+1,1/(1+r^(1/α))];

phinlA[d_,α_,ν_,0,l_,r_] := With[{μ = α(d-2+2*l)}, μ*r^(-l-d+2)*Beta[r^(1/α)/(1+r^(1/α)),μ,ν]];

rhonlA[d_,α_,ν_,0,l_,r_] := With[{μ = α*(d-2+2*l), z = r^(1/(2α))}, r^(l-2+1/α)/(1+z^2)^(μ+ν+1)*((2μ+2ν-1)*(μ+ν))];

phinlA[d_,α_,ν_,n_,l_,r_] := With[{μ = α*(d-2+2*l), z = r^(1/(2α))}, phinlA[d,α,ν,n-1,l,r] - 2*(n-1)!/Pochhammer[μ+1,n-1] * r^l/(1+z^2)^(μ+ν) * JacobiP[n-1,μ+2ν-1,μ,(z^2-1)/(z^2+1)]];

rhonlA[d_,α_,ν_,n_,l_,r_] := With[{μ = α*(d-2+2*l), z = r^(1/(2α))}, r^(l-2+1/α)/(1+z^2)^(μ+ν+1)*((n+2μ+2ν-1)*(n+μ+ν)*JacobiP[n,μ+2ν-1,μ,(z^2-1)/(z^2+1)] - (n+μ+2ν-1)*(n+μ+ν-1)*JacobiP[n-1,μ+2ν-1,μ,(z^2-1)/(z^2+1)])];

NnlA[2,α_,ν_,0,0] := (2*ν-1)^2*α^3*(PolyGamma[0,2*ν] - PolyGamma[0,ν]);

NnlA[d_,α_,ν_,n_,l_] := With[{μ = α*(d-2+2*l)}, α*Gamma[n+μ+2ν]*Gamma[μ+1]/Gamma[n+2μ+2ν-1]];

KnlA[2,α_,ν_,0,0] = 1/(4*Pi);

KnlA[d_,α_,ν_,n_,l_] := With[{μ = α*(d-2+2*l)}, -(n!)*Gamma[μ+1]/(4π*α^2*(2n+2μ+2ν-1)*Gamma[n+μ])];

(* the 'd' parameter in the previous expressions gives the number of
spatial dimensions, so we're only actually interested in d=3: *)

rhonl[α_,ν_,n_,l_,r_] := KnlA[3,α,ν,n,l]*rhonlA[3,α,ν,n,l,r];

phinl[α_,ν_,n_,l_,r_] := phinlA[3,α,ν,n,l,r];

Knl[α_,ν_,n_,l_] := KnlA[3,α,ν,n,l];

ωl[α_, ν_, l_, s_] := Re[2*MellinTransform[rhonl[α,ν,0,l,r], r, 5/2 + I*s] * MellinTransform[rhonl[α,ν,0,l,r], r, 5/2 - I*s]/((l+1/2)^2 + s^2)]; (* polynomial weight function Eq (4.7) *)

PnlA[α_,ν_,n_,l_,s_] := With[{μ = α(1+2*l)}, I^n*(n!) * Gamma[μ]*Gamma[μ+ν]/(Gamma[n+μ]*Gamma[μ+ν+n])*pn[n,μ/2,μ/2+ν,α*s]]; (* the index-raising polynomial (Eq. (4.8)) *)

Pζrho[α_,ν_,n_,l_,r_] := Sum[Coefficient[PnlA[α,ν,n,l,s], s, j]*Nest[ζ[x],rhonl[α,ν,0,l,x],j], {j,0,n}] /. x -> r;

(* Now test with e.g. 
   FullSimplify[rhonl[α,ν,n,l,r]/Pζrho[α,ν,n,l,r]]
   for various values of α, ν, n, l
*)



