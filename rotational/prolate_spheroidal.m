Import["rotational.m"]; Import["operators.m"];

(* coord definition *)

q1 = η; q2 = θ; q3 = φ;

xs = b*Sinh[η]*Sin[θ]*Cos[φ]; ys = b*Sinh[η]*Sin[θ]*Sin[φ]; zs = b*Cosh[η]*Cos[θ];

IJac = Simplify[Inverse[JacMat[xs,ys,zs,η,θ,φ]]];

(* fn, hn, R from Morse & Feshbach *)

h1[η_,θ_] := b*Sqrt[Cosh[η]^2 - Cos[θ]^2];

h2[η_,θ_] := -b*Sqrt[Cosh[η]^2 - Cos[θ]^2];

h3[η_,θ_] := -b*Sinh[η]*Sin[θ];

f1[η_] := b*Sinh[η]; f2[θ_] := -Sin[θ];

R[η_,θ_] := 1;

(* Eigenfunctions: *)

Φalm[η_,θ_,φ_] := -4*Pi/(a^2 + (l+1/2)^2)*Exp[I*m*φ]*LegendreP[I*a-1/2,m,Cosh[η]]*LegendreP[l,m,Cos[θ]];

Ψalm[η_,θ_,φ_] := (g[η,θ])^(-1)*Exp[I*m*φ]*LegendreP[I*a-1/2,m,Cosh[η]]*LegendreP[l,m,Cos[θ]];

(* The operators as defined via symmetries *)

s1[η_,θ_,φ_] := ( dds[dds[#]] - b^2*P3s[P3s[#]] - b^2*Sin[θ]^2*Laps[#] + 1/4*# )&;

s2[η_,θ_,φ_] := ( -dds[dds[#]] + b^2*P3s[P3s[#]] - b^2*Sinh[η]^2*Laps[#] - 1/4*# )&;

S1[η_,θ_,φ_] := ( DDs[DDs[#]] - b^2*P3s[P3s[#]] - b^2*Laps[Sin[θ]^2*#] + 1/4*# )&;

S2[η_,θ_,φ_] := ( -DDs[DDs[#]] + b^2*P3s[P3s[#]] - b^2*Laps[Sinh[η]^2*#] - 1/4*# )&;

(* Check all the operators obey correct relations *)

(* Simplify[(L1[η,θ][f[η,θ,φ]] + L2[η,θ][f[η,θ,φ]] + L3[η,θ,φ][f[η,θ,φ]]) - Laps[f[η,θ,φ]]]

Simplify[(t1[η,θ,φ][f[η,θ,φ]] + t2[η,θ,φ][f[η,θ,φ]]) + g[η,θ]*Laps[f[η,θ,φ]]]

FullSimplify[Comm[t1[η,θ,φ],t2[η,θ,φ]][f[η,θ,φ]]]

FullSimplify[Chop[t1[η,θ,φ][Φalm[η,θ,φ]]/Φalm[η,θ,φ] /. {a -> 4, l -> 3, m -> 2, η -> 3.32, θ -> 1.43, φ -> 0.42}]]

FullSimplify[Chop[t2[η,θ,φ][Φalm[η,θ,φ]]/Φalm[η,θ,φ] /. {a -> 4, l -> 3, m -> 2, η -> 3.32, θ -> 1.43, φ -> 0.42}]]

FullSimplify[Chop[Laps[Φalm[η,θ,φ]]/(4*Pi*Ψalm[η,θ,φ]) /. {a -> 4, l -> 3, m -> 2, η -> 3.32, θ -> 1.234, φ -> 2.22}]] *)

(* alternative (probably numerically unhelpful) eigenfunctions via Whipple's formula *)

Φalm2[η_,θ_,φ_] := -4*Pi/(a^2 + (l+1/2)^2)*1/Sqrt[Sinh[η]]*Exp[I*m*φ]*LegendreP[-m-1/2,-I*a,Coth[η]]*LegendreP[l,m,Cos[θ]];

Ψalm2[η_,θ_,φ_] := (g[η,θ])^(-1)*1/Sqrt[Sinh[η]]*Exp[I*m*φ]*LegendreP[-m-1/2,-I*a,Coth[η]]*LegendreP[l,m,Cos[θ]];


(* Functions for expressing prol. sph. coords in cartesian coords *)

cη2[b_,x_,y_,z_] := 1/2*(1 + (x^2+y^2+z^2)/b^2) + 1/2*Sqrt[1 + 2*(x^2+y^2-z^2)/b^2 + (x^2+y^2+z^2)^2/b^4];

ηxyz[b_,x_,y_,z_] := ArcCosh[Sqrt[cη2[b,x,y,z]]];


