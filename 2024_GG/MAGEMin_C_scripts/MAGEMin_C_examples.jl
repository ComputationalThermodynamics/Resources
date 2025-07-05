# GG2024 - MAGEMin_C - Examples
# Note that: “ig”, igneous; “mp”, metapelite; “mb”, metabasite; “um”, ultramafic

using MAGEMin_C

#= ============================================================== =#
# Example 1 - Single-point minimization with predefined bulk rock
#= ============================================================== =#

data        =   Initialize_MAGEMin("ig", verbose=true);
test        =   0         # KLB-1
data        =   use_predefined_bulk_rock(data, test);
P           =   8.0
T           =   800.0
out         =   single_point_minimization(P,T, data);

Finalize_MAGEMin(data)

#= ============================================================== =#
# Example 2 - Single-point minimization with user defined bulk rock
#= ============================================================== =#

data    = Initialize_MAGEMin("ig", verbose=false);
    
P,T     = 10.0, 1100.0
Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
X       = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
sys_in  = "wt"    
out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)

Finalize_MAGEMin(data)

#= ============================================================== =#
# Example 3 - Multi-point minimization
# Note that the computation will be performed in parallel if Julia is started with multiple threads
#= ============================================================== =#

data    = Initialize_MAGEMin("ig", verbose=false);
P       = [10.0, 12.0]
T       = [1100.0, 1000.0]
Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
X1      = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
X2      = [49.43; 14.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 0.0];
X       = [X1,X2]
sys_in  = "wt"    
out     = multi_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)

Finalize_MAGEMin(data)

#= ============================================================== =#
# Example 4 - Multi-point minimization with pre-allocated output structure
# Note that the computation will be performed in parallel if Julia is started with multiple threads
#= ============================================================== =#

data    = Initialize_MAGEMin("ig", verbose=false);
P       = [10.0, 12.0]
T       = [1100.0, 1000.0]
Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
X1      = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
X2      = [49.43; 14.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 0.0];
X       = [X1,X2]
sys_in  = "wt"    

out     = Vector{MAGEMin_C.gmin_struct{Float64, Int64}}(undef, length(X))
out     = deepcopy(multi_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in))

Finalize_MAGEMin(data)