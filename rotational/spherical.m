Import["rotational.m"]; Import["operators.m"];

(* coord definition *)

q1 = r; q2 = θ; q3 = φ;

xs = r*Sin[θ]*Cos[φ]; ys = r*Sin[θ]*Sin[φ]; zs = r*Cos[θ];

IJac = Simplify[Inverse[JacMat[xs,ys,zs,r,θ,φ]]];


(* fn, hn, R from Morse & Feshbach *)

h1[r_,θ_] := 1; h2[r_,θ_] := -r; h3[r_,θ_] := -r*Sin[θ];

f1[r_] := r^2; f2[θ_] := -Sin[θ];

R[r_,θ_] := 1;


(* Test it with:

Simplify[t1[r,θ,φ][f[r,θ,φ]] - (dds[dds[f[r,θ,φ]]] + 1/4*f[r,θ,φ])]

Simplify[t2[r,θ,φ][f[r,θ,φ]] - (JJs[f[r,θ,φ]])]

Simplify[t1[r,θ,φ][f[r,θ,φ]] + t2[r,θ,φ][f[r,θ,φ]] + g[r,θ]*Laps[f[r,θ,φ]]]

   *)
