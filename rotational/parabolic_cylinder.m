Import["rotational.m"]; Import["operators.m"];

(* coord definition *)

q1 = λ; q2 = μ; q3 = z;

xs = 1/2*(λ^2 - μ^2); ys = λ*μ; zs = z;

IJac = Simplify[Inverse[JacMat[xs,ys,zs,λ,μ,z]]];

h1[λ_,μ_] := Sqrt[λ^2 + μ^2];

h2[λ_,μ_] := Sqrt[λ^2 + μ^2];

h3[λ_,μ_] := 1;

f1[λ_] := 1; f2[μ_] := 1; R[λ_,μ_] := 1;

(* Eigenfunctions: *)

ParabolicCylinderU[a_,z_] := ParabolicCylinderD[-a-1/2,z];

Φαβk[λ_,μ_,z_] := -4*Pi/(α^2 + β^2)*ParabolicCylinderU[(α^2)/(2*k),(2*k)^(1/2)*λ]*ParabolicCylinderU[(β^2)/(2*k),(2*k)^(1/2)*μ]*Exp[I*k*z];

Ψαβk[λ_,μ_,z_] := 1/g[λ,μ]*ParabolicCylinderU[(α^2)/(2*k),(2*k)^(1/2)*λ]*ParabolicCylinderU[(β^2)/(2*k),(2*k)^(1/2)*μ]*Exp[I*k*z];


(* Test it with:

Simplify[(L1[λ,μ][f[λ,μ,z]] + L2[λ,μ][f[λ,μ,z]] + L3[λ,μ,z][f[λ,μ,z]]) - Laps[f[λ,μ,z]]]

Simplify[(t1[λ,μ,z][f[λ,μ,z]] + t2[λ,μ,z][f[λ,μ,z]]) + g[λ,μ]*Laps[f[λ,μ,z]]]

Simplify[Comm[t1[λ,μ,z],t2[λ,μ,z]][f[λ,μ,z]], Assumptions -> {λ>0,μ>0}]

FullSimplify[t1[λ,μ,z][Φαβk[λ,μ,z]]/Φαβk[λ,μ,z] /. {k -> 9/2, μ -> 3, α -> 7}]

FullSimplify[t2[λ,μ,z][Φαβk[λ,μ,z]]/Φαβk[λ,μ,z] /. {k -> 1, λ -> 2, β -> 3}]

FullSimplify[Laps[Φαβk[λ,μ,z]]/(4*Pi*Ψαβk[λ,μ,z]) /. {α -> 9, β -> 4, k -> 5/2}]

*)

