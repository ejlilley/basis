JacMat[x_,y_,z_,q1_,q2_,q3_] := {{D[x,q1], D[x,q2], D[x,q3]}, {D[y,q1], D[y,q2], D[y,q3]}, {D[z,q1], D[z,q2], D[z,q3]}};

p1[q1_] := -h1[q1,q2]*h3[q1,q2]/(R[q1,q2]^2*f1[q1]*f2[q2]*h2[q1,q2]);

p2[q2_] := -h2[q1,q2]*h3[q1,q2]/(R[q1,q2]^2*f1[q1]*f2[q2]*h1[q1,q2]);

g[q1_,q2_] := -R[q1,q2]^2*f1[q1]*f2[q2]*h1[q1,q2]*h2[q1,q2]/(h3[q1,q2]);

hf[q1_,q2_] := h1[q1,q2]*h2[q1,q2]*h3[q1,q2]/(R[q1,q2]^2*f1[q1]*f2[q2]*(-1));

M1[q2_] := Simplify[hf[x1,x2]/h1[x1,x2]^2] /. x2 -> q2;

M2[q1_] := Simplify[hf[x1,x2]/h2[x1,x2]^2] /. x1 -> q1;

t3[q3_] := ( -( D[#,{q3,2}] ) ) &;

(* the "x0" thing is a hack to get correct behaviour when the expression doesn't depend on one of the variables *)

v1[q1_,q2_] := (Select[x0 + Expand[FullSimplify[g[x1,x2]/h3[x1,x2]^2]],FreeQ[x2]]) /. {x1 -> q1, x2 -> q2, x0 -> 0};

v2[q1_,q2_] := (Select[x0 + Expand[FullSimplify[g[x1,x2]/h3[x1,x2]^2]],FreeQ[x1]]) /. {x1 -> q1, x2 -> q2, x0 -> 0};

(* v1[q1_] := Select[x0 + FullSimplify[g[x1,x2]/h3[x1,x2]^2],FreeQ[x2]] /. {x1 -> q1, x0 -> 0};

v2[q2_] := Select[x0 + FullSimplify[g[x1,x2]/h3[x1,x2]^2],FreeQ[x1]] /. {x2 -> q2, x0 -> 0}; *)

u[q1_,q2_] := (FullSimplify[-1/2*Integrate[FullSimplify[p2[x2]*(D[Log[R[x1,x2]^2*f1[x1]],x1]*D[D[Log[R[x1,x2]^2],x1],x2] + D[D[Log[R[x1,x2]^2],x2],{x1,2}])],x2]]) /. {x1 -> q1, x2 -> q2};

(* Alternatively: *)

u2[q1_,q2_] := (FullSimplify[-1/2*Integrate[FullSimplify[p1[x2]*(D[Log[R[x1,x2]^2*f2[x2]],x2]*D[D[Log[R[x1,x2]^2],x2],x1] + D[D[Log[R[x1,x2]^2],x1],{x2,2}])],x1]]) /. {x1 -> q1, x2 -> q2};



(* The actual operators: *)

(* t1[q1_,q2_,q3_] := ( 1/(M2[q1]) * ( 1/(R[q1,q2]^2*f1[q1])*D[R[q1,q2]^2*f1[q1],q1]*D[#,q1] + D[#,{q1,2}] ) + v1[q1,q2]*D[#, {q3,2}] ) &;

t2[q1_,q2_,q3_] := ( 1/(M1[q2]) * ( 1/(R[q1,q2]^2*f2[q2])*D[R[q1,q2]^2*f2[q2],q2]*D[#,q2] + D[#,{q2,2}] ) + v2[q1,q2]*D[#, {q3,2}] ) &; *)

t1[q1_,q2_,q3_] := ( -(g[q1,q2]/(h1[q1,q2])^2 * ( 1/(R[q1,q2]^2*f1[q1])*D[R[q1,q2]^2*f1[q1],q1]*D[#,q1] + D[#,{q1,2}] ) + v1[q1,q2]*D[#, {q3,2}] - u[q1,q2]*#) ) &;

t2[q1_,q2_,q3_] := ( -(g[q1,q2]/(h2[q1,q2])^2 * ( 1/(R[q1,q2]^2*f2[q2])*D[R[q1,q2]^2*f2[q2],q2]*D[#,q2] + D[#,{q2,2}] ) + v2[q1,q2]*D[#, {q3,2}] + u[q1,q2]*#) ) &;

L1[q1_,q2_] := ( 1/(h1[q1,q2])^2 * ( 1/(R[q1,q2]^2*f1[q1])*D[R[q1,q2]^2*f1[q1],q1]*D[#,q1] + D[#,{q1,2}] ) ) &;

L2[q1_,q2_] := ( 1/(h2[q1,q2])^2 * ( 1/(R[q1,q2]^2*f2[q2])*D[R[q1,q2]^2*f2[q2],q2]*D[#,q2] + D[#,{q2,2}] ) ) &;

L3[q1_,q2_,q3_] := ( ( 1/(h3[q1,q2])^2 * D[#,{q3,2}] ) ) &;
