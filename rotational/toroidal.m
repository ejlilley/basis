Import["rotational.m"]; Import["operators.m"];

(* coord definition *)

q1 = μ; q2 = θ; q3 = φ;

xs = b*Sinh[μ]*Cos[φ]/(Cosh[μ] - Cos[θ]); ys = b*Sinh[μ]*Sin[φ]/(Cosh[μ] - Cos[θ]); zs = b*Sin[θ]/(Cosh[μ] - Cos[θ]);

IJac = Simplify[Inverse[JacMat[xs,ys,zs,μ,θ,φ]]];

(* fn, hn, R from Morse & Feshbach *)

h1[μ_,θ_] := b/(Cosh[μ] - Cos[θ]);

h2[μ_,θ_] := -b/(Cosh[μ] - Cos[θ]);

h3[μ_,θ_] := -b*Sinh[μ]/(Cosh[μ] - Cos[θ]);

f1[μ_] := Sinh[μ]; f2[θ_] := -1; R[μ_,θ_] := Sqrt[b/(Cosh[μ] - Cos[θ])];

(* Eigenfunctions: *)

Φslm[μ_,θ_,φ_] := -4*Pi/R[μ,θ]*LegendreP[I*s-1/2,m,Cosh[μ]]*Exp[I*m*φ+I*l*θ]/(s^2 + l^2);

Ψslm[μ_,θ_,φ_] := 1/R[μ,θ]^5*LegendreP[I*s-1/2,m,Cosh[μ]]*Exp[I*m*φ+I*l*θ];


(* Alternative forms via Whipple formula *)

Φslm2[μ_,θ_,φ_] := -4*Pi/(Sqrt[Sinh[μ]]*R[μ,θ])*LegendreP[-m-1/2,-I*s,Coth[μ]]*Exp[I*m*φ+I*l*θ]/(s^2 + l^2);

Ψslm2[μ_,θ_,φ_] := 1/(Sqrt[Sinh[μ]]*R[μ,θ]^5)*LegendreP[-m-1/2,-I*s,Coth[μ]]*Exp[I*m*φ+I*l*θ];


(* Expected oeprators: *)

s1 := ( -(1/b^2*Sq[k3s][#] + 1/2*Symm[k3s,P3s][#] + b^2/4*Sq[P3s][#]) - g[μ,θ]*Laps[#] + 1/4*#)&;

s2s := ( 1/b*k3s[#] + b/2*P3s[#])&;

s2 := ( s2s[s2s[#]] - 1/4*#)&;


(* Test it with:

Simplify[(L1[μ,θ][f[μ,θ,φ]] + L2[μ,θ][f[μ,θ,φ]] + L3[μ,θ,φ][f[μ,θ,φ]]) - Laps[f[μ,θ,φ]]]

Simplify[(t1[μ,θ,φ][f[μ,θ,φ]] + t2[μ,θ,φ][f[μ,θ,φ]]) + g[μ,θ]*Laps[f[μ,θ,φ]]]

Simplify[Comm[t1[μ,θ,φ],t2[μ,θ,φ]][f[μ,θ,φ]], Assumptions -> {μ>0,θ>0}]

Chop[FullSimplify[t1[μ,θ,φ][Φslm[μ,θ,φ]]/(Φslm[μ,θ,φ]) /. {m -> 2, s -> 4, l -> 3, μ -> 1.2, θ -> 0.123}] - 1/4]

Chop[FullSimplify[t2[μ,θ,φ][Φslm[μ,θ,φ]]/(Φslm[μ,θ,φ]) /. {m -> 2, s -> 4, l -> 3, μ -> 1.2, θ -> 0.123}] + 1/4]

Chop[N[FullSimplify[Laps[Φslm[μ,θ,φ]]/(4*Pi*Ψslm[μ,θ,φ]) /. {θ -> 0.1, φ -> 0.2, μ -> 1.3}] /. {s -> 14, l -> 4, m -> 2}]]

Simplify[(t1[μ,θ,φ][f[μ,θ,φ]]) - s1[f[μ,θ,φ]]]

Simplify[(t2[μ,θ,φ][f[μ,θ,φ]]) - s2[f[μ,θ,φ]]]

   *)
