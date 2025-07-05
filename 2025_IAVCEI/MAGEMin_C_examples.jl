using MAGEMin_C
data    = Initialize_MAGEMin("ig", verbose=false; solver=0);
P,T     = 10.0, 1100.0
Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
X       = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
sys_in  = "wt"
out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)



#=

typeof(out) = MAGEMin_C.gmin_struct{Float64, Int64}
out_XY = Vector{MAGEMin_C.gmin_struct{Float64, Int64}}(undef,2)

=#