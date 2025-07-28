Import["rotational.m"]; Import["operators.m"];

q1 = r; q2 = φ; q3 = z;

xs = r*Cos[φ]; ys = r*Sin[φ]; zs = z;

IJac = Simplify[Inverse[JacMat[xs,ys,zs,r,φ,z]]];

h1[r_,φ_] := 1; h2[r_,φ_] := -r; h3[r_,φ_] := -1; (* note Morse & Feshbach has h3 = 1 *)

f1[r_] := r; f2[φ_] := -1; R[r_,φ_] := 1;

(* Eigenfunctions: *)

Φskm[r_,φ_,z_] := -4*Pi/(s^2 + m^2)*Exp[I*m*φ + I*k*z]*BesselK[I*s,k*r];

Ψskm[r_,φ_,z_] := 1/r^2*Exp[I*m*φ + I*k*z]*BesselK[I*s,k*r];

(* Test it with:

Simplify[(L1[r,φ][f[r,φ,z]] + L2[r,φ][f[r,φ,z]] + L3[r,φ,z][f[r,φ,z]]) - Laps[f[r,φ,z]]]

Simplify[(t1[r,φ,z][f[r,φ,z]] + t2[r,φ,z][f[r,φ,z]]) + g[r,φ]*Laps[f[r,φ,z]]]

FullSimplify[Comm[t1[r,φ,z],t2[r,φ,z]][f[r,φ,z]]]

FullSimplify[(t1[r,φ,z][Φskm[r,φ,z]])/Φskm[r,φ,z]]

FullSimplify[(Laps[Φskm[r,φ,z]])/(4*Pi*Ψskm[r,φ,z])]

FullSimplify[(t2[r,φ,z][f[r,φ,z]])]

   *)
