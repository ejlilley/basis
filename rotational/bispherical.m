Import["rotational.m"]; Import["operators.m"];

(* coord definition *)

q1 = μ; q2 = θ; q3 = φ;

xs = b*Sin[θ]*Cos[φ]/(Cosh[μ] - Cos[θ]); ys = b*Sin[θ]*Sin[φ]/(Cosh[μ] - Cos[θ]); zs = b*Sinh[μ]/(Cosh[μ] - Cos[θ]);

IJac = Simplify[Inverse[JacMat[xs,ys,zs,μ,θ,φ]]];

(* fn, hn, R from Morse & Feshbach *)

(* g[q1_,q2_] := -R[q1,q2]^2*f1[q1]*f2[q2]*h1[q1,q2]*h2[q1,q2]/(h3[q1,q2]); *)

h1[μ_,θ_] := b/(Cosh[μ] - Cos[θ]);

h2[μ_,θ_] := -b/(Cosh[μ] - Cos[θ]);

h3[μ_,θ_] := -b*Sin[θ]/(Cosh[μ] - Cos[θ]);

f1[μ_] := 1; f2[θ_] := -Sin[θ]; R[μ_,θ_] := Sqrt[b/(Cosh[μ] - Cos[θ])];

(* Eigenfunctions: *)

Φslm[μ_,θ_,φ_] := 4*Pi*b^2*(R[μ,θ])^(-1)/((s)^2 + (l+1/2)^2)*Exp[I*s*μ]*LegendreP[l,m,Cos[θ]]*Exp[I*m*φ];

Ψslm[μ_,θ_,φ_] := (R[μ,θ])^(-5)*Exp[I*s*μ]*LegendreP[l,m,Cos[θ]]*Exp[I*m*φ];

(* Expected operators: *)

s1s[μ_,θ_,φ_] := ( 1/b*k3s[#] - b/2*P3s[#])&;

s1[μ_,θ_,φ_] := ( Sq[s1s[μ,θ,φ]][#] )&;

s2[μ_,θ_,φ_] := ( -s1[μ,θ,φ][#] - 1/4*# - b^2/(Cosh[μ] - Cos[θ])^2*Laps[#] )&;

(* Test it with: *)

(* Simplify[(L1[μ,θ][f[μ,θ,φ]] + L2[μ,θ][f[μ,θ,φ]] + L3[μ,θ,φ][f[μ,θ,φ]]) - Laps[f[μ,θ,φ]]]

Simplify[(t1[μ,θ,φ][f[μ,θ,φ]] + t2[μ,θ,φ][f[μ,θ,φ]]) + g[μ,θ]*Laps[f[μ,θ,φ]]]

Simplify[Comm[t1[μ,θ,φ],t2[μ,θ,φ]][f[μ,θ,φ]], Assumptions -> {μ>0,θ>0}]

Simplify[(s1[μ,θ,φ][f[μ,θ,φ]] + s2[μ,θ,φ][f[μ,θ,φ]] + 1/4*f[μ,θ,φ]) + g[μ,θ]*Laps[f[μ,θ,φ]]]

FullSimplify[t2[μ,θ,φ][Φslm[μ,θ,φ]]/Φslm[μ,θ,φ] /. {s -> 1, l -> 4, m -> 1}, Assumptions -> {μ>0,θ>0}] - 1/4

FullSimplify[t1[μ,θ,φ][Φslm[μ,θ,φ]]/Φslm[μ,θ,φ], Assumptions -> {μ>0,θ>0}]

FullSimplify[t1[μ,θ,φ][f[μ,θ,φ]] - s1[μ,θ,φ][f[μ,θ,φ]]]

FullSimplify[t2[μ,θ,φ][f[μ,θ,φ]] - s2[μ,θ,φ][f[μ,θ,φ]] - 1/4*f[μ,θ,φ]]

   *)
