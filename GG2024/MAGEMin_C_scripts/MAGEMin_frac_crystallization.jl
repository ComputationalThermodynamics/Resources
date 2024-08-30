# GG2024 - MAGEMin_C - Fractional crystallization
# Note that: “ig”, igneous; “mp”, metapelite; “mb”, metabasite; “um”, ultramafic

using MAGEMin_C
using Plots

dtb     = "ig"
data    = Initialize_MAGEMin(dtb);


Xoxides = ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "O", "Cr2O3", "H2O"]
X       = [50.0, 8.7, 11.7, 12.14, 7.7, 0.2, 2.5, 1.0, 0.5, 0.01, 10.0]
sys_in  = "mol";


P       = 5.0;                                                              # Pressure in kbar                                     
n_steps = 64;                                                               # Number of steps in the Fractional Crystallization path               
Ts      = 1200.0;                                                           # Starting temperature in the Fractional Crystallization path           
Te      = 600.0;                                                            # Ending temperature in the Fractional Crystallization path     
T       = Array(range(Ts, stop=Te, length=n_steps))                         # Defines temperature values for the Fractional Crystallization path
out     = Vector{MAGEMin_C.gmin_struct{Float64, Int64}}(undef, n_steps)     # Vector to store the output of the single_point_minimization function

#= Few comments 
    Here deepcopy is used to copy the output of the single_point_minimization and not only the reference to the output.
    This is important because the output of the single_point_minimization function is a mutable struct and if we do not use deepcopy,
    the output of the single_point_minimization function will be overwritten in the next iteration of the loop.
=#
for i in 1:n_steps
    out[i] = deepcopy(single_point_minimization(P, T[i], data, X=X, Xoxides=Xoxides, sys_in=sys_in))
    if "liq" in out[i].ph                                                   # If the liquid phase is present, the bulk composition is updated
        X = deepcopy(out[i].bulk_M)
    end
    # Note that the bulk composition is updated only if the liquid phase is present
    # This is because the bulk composition is the composition of the liquid phase at the previous step
    # If the liquid phase is not present, the bulk composition is not updated and the composition of the liquid phase at the previous step is used
end

Finalize_MAGEMin(data)

#=
    In the following section we extract the melt fraction, total melt fraction, SiO2 in the melt, melt density for all steps
=#
frac_M      = [out[i].frac_M for i in 1:n_steps];                           # Melt fraction for all steps
frac_M_tot  = accumulate(*, frac_M)

SiO2_id     = findfirst(out[1].oxides .== "SiO2")                           # Index of SiO2 in the oxides array   
dry_id      = findall(out[1].oxides .!= "H2O")                              # Indices of all oxides except H2O
SiO2_M_dry  = [ (out[i].bulk_M[SiO2_id] / sum(out[i].bulk_M[dry_id])*100.0) for i in 1:n_steps];                             # SiO2 in the melt for all steps
rho_M       = [ (out[i].rho_M) for i in 1:n_steps];                         # melt density for all steps

rho_M[rho_M .== 0.0] .= NaN;                                                # Replace 0.0 values with NaN


#=
    Ploting the results using Plots
=#
p1          = plot(T,frac_M,     xflip=true, xlabel="Temperature (°C)", marker = :circle, markersize = 2, lw=2, ylabel="Melt fraction (mol)",         legend=false)
p2          = plot(T,frac_M_tot, xflip=true, xlabel="Temperature (°C)", marker = :circle, markersize = 2, lw=2, ylabel="Total melt fraction (mol)",   legend=false)
p3          = plot(T,rho_M,      xflip=true, xlabel="Temperature (°C)", marker = :circle, markersize = 2, lw=2, ylabel="Melt density (kg/m³)",        legend=false)
p4          = plot(T,SiO2_M_dry, xflip=true, xlabel="Temperature (°C)", marker = :circle, markersize = 2, lw=2, ylabel="SiO₂ melt anhydrous (mol%)",  legend=false)

fig = plot(p1, p2, p3, p4, layout=(2, 2), size=(800, 600))
savefig(fig,"frac_crystallization.png")