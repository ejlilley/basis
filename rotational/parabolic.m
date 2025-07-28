Import["rotational.m"]; Import["operators.m"];

q1 = λ; q2 = μ; q3 = φ;

xs = λ*μ*Cos[φ]; ys = λ*μ*Sin[φ]; zs = 1/2*(λ^2 - μ^2);

IJac = Simplify[Inverse[JacMat[xs,ys,zs,λ,μ,φ]]];

h1[λ_,μ_] := Sqrt[λ^2 + μ^2]; h2[λ_,μ_] := Sqrt[λ^2 + μ^2]; h3[λ_,μ_] := -λ*μ;

f1[λ_] := λ; f2[μ_] := μ; R[λ_,μ_] := 1;

(* Eigenfunctions: *)

Φkqm[λ_,μ_,φ_] := -4*Pi/(k^2 + q^2)*BesselJ[m,k*λ]*BesselJ[m,q*μ]*Exp[I*m*φ];

Ψkqm[λ_,μ_,φ_] := 1/g[λ,μ]*BesselJ[m,k*λ]*BesselJ[m,q*μ]*Exp[I*m*φ];

(* Operators in terms of symmetries: *)

s1 := ( 1/2*(Symm[P1s,J2s][#] - Symm[P2s,J1s][#] + Symm[P3s,dds][#] - (λ^2+μ^2)*Laps[#]) ) &;

s2 := ( -1/2*(Symm[P1s,J2s][#] - Symm[P2s,J1s][#] + Symm[P3s,dds][#] + (λ^2+μ^2)*Laps[#]) ) &;

(* Test it with:

Simplify[(L1[λ,μ][f[λ,μ,φ]] + L2[λ,μ][f[λ,μ,φ]] + L3[λ,μ,φ][f[λ,μ,φ]]) - Laps[f[λ,μ,φ]]]

Simplify[(t1[λ,μ,φ][f[λ,μ,φ]] + t2[λ,μ,φ][f[λ,μ,φ]]) + g[λ,μ]*Laps[f[λ,μ,φ]]]

Simplify[Comm[t1[λ,μ,φ],t2[λ,μ,φ]][f[λ,μ,φ]], Assumptions -> {λ>0,μ>0}]

FullSimplify[Laps[Φkqm[λ,μ,φ]]/(4*Pi*Ψkqm[λ,μ,φ])]

FullSimplify[t1[λ,μ,φ][Φkqm[λ,μ,φ]]/Φkqm[λ,μ,φ]]

FullSimplify[t2[λ,μ,φ][Φkqm[λ,μ,φ]]/Φkqm[λ,μ,φ]]

FullSimplify[t1[λ,μ,φ][f[λ,μ,φ]] - s1[f[λ,μ,φ]]]

FullSimplify[t2[λ,μ,φ][f[λ,μ,φ]] - s2[f[λ,μ,φ]]]

*)
