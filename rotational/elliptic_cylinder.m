Import["rotational.m"]; Import["operators.m"];

(* coord definition *)

q1 = μ; q2 = φ; q3 = z;

xs = b*Cosh[μ]*Cos[φ]; ys = b*Sinh[μ]*Sin[φ]; zs = z;

IJac = Simplify[Inverse[JacMat[xs,ys,zs,μ,φ,z]]];

h1[μ_,φ_] := b*Sqrt[Cosh[μ]^2 - Cos[φ]^2];

h2[μ_,φ_] := -b*Sqrt[Cosh[μ]^2 - Cos[φ]^2];

h3[μ_,φ_] := 1;

f1[μ_] := 1; f2[φ_] := -1; R[μ_,φ_] := 1;

(* Eigenfunctions: *)

q = (b*k/2)^2; me[ν_,z_,q_] := MathieuC[ν,q,z] + I*MathieuS[ν,q,z];

Φamk[μ_,φ_,z_] := -4*Pi/(a^2 + m^2)*Exp[I*k*z]*me[2*q + a^2, I*μ - Pi/2, q]*me[2*q - m^2, φ - Pi/2, q];

Ψamk[μ_,φ_,z_] := 1/g[μ,φ]*Exp[I*k*z]*me[2*q + a^2, I*μ - Pi/2, q]*me[2*q - m^2, φ - Pi/2, q];


(* Test with:

Simplify[(L1[μ,φ][f[μ,φ,z]] + L2[μ,φ][f[μ,φ,z]] + L3[μ,φ,z][f[μ,φ,z]]) - Laps[f[μ,φ,z]]]

Simplify[(t1[μ,φ,z][f[μ,φ,z]] + t2[μ,φ,z][f[μ,φ,z]]) + g[μ,φ]*Laps[f[μ,φ,z]]]

FullSimplify[Comm[t1[μ,φ,z],t2[μ,φ,z]][f[μ,φ,z]]]

FunctionExpand[FullSimplify[t2[μ,φ,z][Φamk[μ,φ,z]]/Φamk[μ,φ,z]]]

Simplify[Laps[Φamk[μ,φ,z]]/(4*Pi*Ψamk[μ,φ,z])]

   *)



