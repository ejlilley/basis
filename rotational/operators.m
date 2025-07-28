K1 := (I*(5/2*x*# + 1/2*(x^2 - y^2 - z^2)*D[#,x] + x*z*D[#,z] + x*y*D[#,y] ))&;

K2 := (I*(5/2*y*# + 1/2*(y^2 - z^2 - x^2)*D[#,y] + y*x*D[#,x] + y*z*D[#,z] ))&;

K3 := (I*(5/2*z*# + 1/2*(z^2 - x^2 - y^2)*D[#,z] + z*y*D[#,y] + z*x*D[#,x] ))&;

k1 := (K1[#] - 2*I*x*#)&;

k2 := (K2[#] - 2*I*y*#)&;

k3 := (K3[#] - 2*I*z*#)&;

DD := (I*(x*D[#,x] + y*D[#,y] + z*D[#,z] + 5/2*#))&;

dd := (I*(x*D[#,x] + y*D[#,y] + z*D[#,z] + 1/2*#))&;

Lap := (D[#,{x,2}] + D[#,{y,2}] + D[#,{z,2}])&;

J1 := (I*(-y*D[#,z] + z*D[#,y]))&; J2 := (I*(-z*D[#,x] + x*D[#,z]))&; J3 := (I*(-x*D[#,y] + y*D[#,x]))&;

P1 := I*D[#,x]&; P2 := I*D[#,y]&; P3 := I*D[#,z]&;

JJ := (J1[J1[#]] + J2[J2[#]] + J3[J3[#]])&;

(* brackets of operators *)

Comm[s_,t_] := (s[t[#]] - t[s[#]])&;

Symm[s_,t_] := (s[t[#]] + t[s[#]])&;

Add[s_,t_] := (s[#] + t[#])&;

Sub[s_,t_] := (s[#] - t[#])&;

Mul[s_,t_] := (s[t[#]])&;

Sq[s_] := (s[s[#]])&;

(* Basis of conformal Killing vectors *)

P1v = {1,0,0}; P2v = {0,1,0}; P3v = {0,0,1};

J1v = {0,z,-y}; J2v = {-z,0,x}; J3v = {y,-x,0};

k1v = {1/2*(x^2-y^2-z^2), x*y, x*z}; k2v = {y*x, 1/2*(y^2-x^2-z^2), y*z}; k3v = {z*x, z*y, 1/2*(z^2-y^2-x^2)};

dv = {x,y,z};

vbasis = {P1v, P2v, P3v, J1v, J2v, J3v, k1v, k2v, k3v, dv};

vbasisnames = {"P1", "P2", "P3", "J1", "J2", "J3", "k1", "k2", "k3", "d"};


(* in different coords defined implicitly via (xs, ys, zs) *)

P1s = (I*(IJac[[1]][[1]]*D[#,q1] + IJac[[2]][[1]]*D[#,q2] + IJac[[3]][[1]]*D[#,q3]))&;

P2s = (I*(IJac[[1]][[2]]*D[#,q1] + IJac[[2]][[2]]*D[#,q2] + IJac[[3]][[2]]*D[#,q3]))&;

P3s = (I*(IJac[[1]][[3]]*D[#,q1] + IJac[[2]][[3]]*D[#,q2] + IJac[[3]][[3]]*D[#,q3]))&;

J1s = ((-ys*P3s[#] + zs*P2s[#]))&; J2s = ((-zs*P1s[#] + xs*P3s[#]))&; J3s = ((-xs*P2s[#] + ys*P1s[#]))&;

JJs = (J1s[J1s[#]] + J2s[J2s[#]] + J3s[J3s[#]])&;

DDs = (xs*P1s[#] + ys*P2s[#] + zs*P3s[#] + I*5/2*#)&;

dds = (xs*P1s[#] + ys*P2s[#] + zs*P3s[#] + I*1/2*#)&;

K1s = (I*5/2*xs*# + 1/2*(xs^2 - ys^2 - zs^2)*P1s[#] + xs*zs*P3s[#] + xs*ys*P2s[#] )&;

K2s = (I*5/2*ys*# + 1/2*(ys^2 - zs^2 - xs^2)*P2s[#] + ys*xs*P1s[#] + ys*zs*P3s[#] )&;

K3s = (I*5/2*zs*# + 1/2*(zs^2 - xs^2 - ys^2)*P3s[#] + zs*ys*P2s[#] + zs*xs*P1s[#] )&;

k1s = (I*1/2*xs*# + 1/2*(xs^2 - ys^2 - zs^2)*P1s[#] + xs*zs*P3s[#] + xs*ys*P2s[#] )&;

k2s = (I*1/2*ys*# + 1/2*(ys^2 - zs^2 - xs^2)*P2s[#] + ys*xs*P1s[#] + ys*zs*P3s[#] )&;

k3s = (I*1/2*zs*# + 1/2*(zs^2 - xs^2 - ys^2)*P3s[#] + zs*ys*P2s[#] + zs*xs*P1s[#] )&;

Laps = ( (-Sq[P1s][#] - Sq[P2s][#] - Sq[P3s][#]) & );
