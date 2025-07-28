Import["rotational.m"]; Import["operators.m"];

q1 = μ; q2 = θ ; q3 = φ;

xs = b*Cosh[μ]*Sin[θ]*Cos[φ]; ys = b*Cosh[μ]*Sin[θ]*Sin[φ]; zs = b*Sinh[μ]*Cos[θ];

IJac = Simplify[Inverse[JacMat[xs,ys,zs,μ,θ,φ]]];

h1[μ_,θ_] := b*Sqrt[Sinh[μ]^2 + Cos[θ]^2];

h2[μ_,θ_] := -b*Sqrt[Sinh[μ]^2 + Cos[θ]^2];

h3[μ_,θ_] := -b*Sin[θ]*Cosh[μ];

f1[μ_] := b*Cosh[μ]; f2[θ_] := -Sin[θ]; R[μ_,θ_] := 1;

(* Eigenfunctions: *)

Φalm[μ_,θ_,φ_] := -4*Pi/(a^2 + (l+1/2)^2)*Exp[I*m*φ]*LegendreP[I*a-1/2,m,I*Sinh[μ]]*LegendreP[l,m,Cos[θ]];

Ψalm[η_,θ_,φ_] := (g[η,θ])^(-1)*Exp[I*m*φ]*LegendreP[I*a-1/2,m,I*Sinh[η]]*LegendreP[l,m,Cos[θ]];


(* Expected operators: *)

s1[μ_,θ_,φ_] := ( dds[dds[#]] + b^2*P3s[P3s[#]] + b^2*Sin[θ]^2*Laps[#] + 1/4*# )&;

s2[μ_,θ_,φ_] := ( -dds[dds[#]] - b^2*P3s[P3s[#]] - b^2*Cosh[μ]^2*Laps[#] - 1/4*# )&;

(* Test it with:

Simplify[(L1[μ,θ][f[μ,θ,φ]] + L2[μ,θ][f[μ,θ,φ]] + L3[μ,θ,φ][f[μ,θ,φ]]) - Laps[f[μ,θ,φ]]]

Simplify[(t1[μ,θ,φ][f[μ,θ,φ]] + t2[μ,θ,φ][f[μ,θ,φ]]) + g[μ,θ]*Laps[f[μ,θ,φ]]]

FullSimplify[Comm[t1[μ,θ,φ],t2[μ,θ,φ]][f[μ,θ,φ]]]

Chop[FullSimplify[t1[μ,θ,φ][Φalm[μ,θ,φ]]/Φalm[μ,θ,φ] /. {a -> 3, l -> 1, m -> 1, μ -> 1.1}]] - 1/4

Chop[FullSimplify[t2[μ,θ,φ][Φalm[μ,θ,φ]]/Φalm[μ,θ,φ] /. {a -> 3, l -> 3, m -> 1, μ -> 1.1, θ -> 2.3}]]

FullSimplify[Chop[Laps[Φalm[μ,θ,φ]]/(4*Pi*Ψalm[μ,θ,φ]) /. {a -> 4, l -> 3, m -> 2, μ -> 3.32, θ -> 1.234, φ -> 2.22}]]


FullSimplify[t1[μ,θ,φ][f[μ,θ,φ]] - s1[μ,θ,φ][f[μ,θ,φ]]]

FullSimplify[t2[μ,θ,φ][f[μ,θ,φ]] - s2[μ,θ,φ][f[μ,θ,φ]]]

   *)
