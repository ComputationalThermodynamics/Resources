## Example of loop phase equilibrium calculation using MAGEMin_C


First let's first initialize MAGEMin with the metapelite database (White et al ., 2014)

```julia
data        =   Initialize_MAGEMin("mp", verbose=true);
```

Then define the bulk-rock composition (wt fraction) the related oxide list and the system unit

```julia
X           = [0.5922, 0.1813, 0.006, 0.0223, 0.0633, 0.0365, 0.0127, 0.0084, 0.0016, 0.0007, 0.075]
Xoxides     = ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "O", "MnO", "H2O"]
sys_unit    = "wt"
````

> [!NOTE]
> Here we use the water-oversaturated FWorld Average composition for metapelite

Define test pressure and temperature condition for test point:
```julia
P = 10.0
T = 700.0
```

and perform the test calculation:

```julia
out  = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_unit)
```

which should gives:

```
Pressure          : 10.0      [kbar]
Temperature       : 700.0    [Celsius]
     Stable phase | Fraction (mol fraction) 
                g   0.07753 
               mu   0.169 
              liq   0.24053 
               bi   0.09972 
              ilm   0.01863 
                q   0.24848 
               ky   0.04752 
               ru   0.00037 
              H2O   0.09821 
     Stable phase | Fraction (wt fraction) 
                g   0.09343 
               mu   0.19964 
              liq   0.20529 
               bi   0.10706 
              ilm   0.02116 
                q   0.27091 
               ky   0.06987 
               ru   0.00054 
              H2O   0.03211 
     Stable phase | Fraction (vol fraction) 
                g   0.0595 
               mu   0.18 
              liq   0.25537 
               bi   0.09142 
              ilm   0.01096 
                q   0.26397 
               ky   0.04914 
               ru   0.00033 
              H2O   0.0893 
Gibbs free energy : -853.149024  (27 iterations; 12.19 ms)
Oxygen fugacity          : -13.51247569580698
Delta QFM                : 2.4956538482297375
```

Instead, of using a single pressure and temperature conditions let's now keep the pressure fixed and vary the temperature

```julia
n = 50
P = 10.0
T = collect(range(400.0,800.0,n))
```

Here ```range(min,max,n)``` will create a range of value between min and max with n steps. ```collect()``` turns the range into an array that can be indexed.


Create a loop from 1 to n in order to compute a stable phase equilibrium for all the entries of the ```T```array. This can be done as follow:

```julia
for i=1:n
    T_calc = T[i]       # retrieves the temperature from the temperature array we just defined  
    out  = single_point_minimization(P, T_calc, data, X=X, Xoxides=Xoxides, sys_in=sys_unit)
end
```

Although this performs the set of 50 phase equilibrium calculations as intended, the script is not yet very useful as the results of the calculation are not stored.

To store the results of the stable phase equilibrium calculations you can declare an array of MAGEMin_C output structure such as:

```julia
out = Vector{MAGEMin_C.gmin_struct{Float64, Int64}}(undef,n)
for i=1:n
    T_calc  = T[i]       # retrieves the temperature from the temperature array we just defined  
    out[i]  = single_point_minimization(P, T_calc, data, X=X, Xoxides=Xoxides, sys_in=sys_unit)
end
```

> [!NOTE]
> The first line of previous snippet create a Vector of MAGEMin_C structures, of size ```n```and with undefined content. In this case we are interested in 1D array, but in a similar manner you could create a Matrix e.g., ```out = Matrix{MAGEMin_C.gmin_struct{Float64, Int64}}(undef,n,n)```

Once the calculation are peformed you can access all the informations:

```
out[2].    #then hit double `tabulation` key twice
```

The latter will display all the entries of point ```2```. If you want to retrieve the melt volume fraction you can do so by accessing ```out[2].frac_M_vol```. In the case of point 2, we get:

```
julia> out[2].frac_M_vol
0.0
```

Now the question is how to gather in a efficient `juliaesk`manner the volume fractions for all computed points? One relatively quick and efficient way is as follow:

```julia
frac_M_vol = [out[i].frac_M_vol for i=1:n]
```

which should yield in the terminal:
```
julia> frac_M_vol = [out[i].frac_M_vol for i=1:n]
50-element Vector{Float64}:
 0.0
 0.0
 0.0
 0.0
 0.0
 ⋮
 0.66483750887145
 0.6823915840521158
 0.7005692654225463
 0.7193913179922915
 0.7388665024379735
```


To retrieve the vol fraction of a solid phase, the process is slightly more complex as we have to check in the phase is present in the mineral assemblage first. Let's first allocate an array for the volume fraction of quartz (`q`):

```julia
q_vol = zeros(Float64,n)
```

Let us now loop through the output structure to look for quartz in the stable mineral assemblage, and in the event quartz is stable, retrieve its volume fraction:

```julia
for i=1:n
    if "q" in out[i].ph     #here we check if out[i].ph contains "q"
        id_q = findfirst(out[i].ph .== "q")
        q_vol[i] = out[i].ph_frac_vol[id_q]
    end
end
```

> [!NOTE]
> We first here check if "q" is in the phase assemblage. Then the command `id_q = findfirst(out[i].ph .== "q")`find the position of "quartz" in the array to be able to retrieve to right value.


### Exercises

1. Using `plots.jl` and the julia introduction, plot the volume fraction of melt and quartz as function of the temperature

2. Add to the tlot the volume fraction evolution of `chl`and `bi`

3. Keep the metapelite composition for now. Given a an average geotherm of 10°C/km for a subduction zone create a depth-dependent pressure-temperature path (using 50 steps as before). The starting conditions are 100°C and 10km depth with a maximum depth of 100km. Use a simple relationship to convert depth to pressure such as:

    ```
    P (GPa) = ρgz / 1e9,
    ```
    where `ρ = 3300.0`, `g = 9.81` and `z` is in meter. Don't forget that pressure in MAGEMin as to be converted to `kbar`


    - Plot the evolution of the volume fraction of liquid and free water "H2O". 
    - What needs to be adjusted to make the modelling more realistic?


<!-- 
```julia

geo = 10.0
z   = collect(range(10.0,100,n))
P   = z .* 1000.0 .* 3300 .* 9.81 ./ 1e8
T   = 100.0 .+ z .* geo

for i=1:n
    T_calc  = T[i]       # retrieves the temperature from the temperature array we just defined  
    P_calc  = P[i]
    out[i]  = single_point_minimization(P_calc, T_calc, data, X=X, Xoxides=Xoxides, sys_in=sys_unit)
end

```
-->

## Example with non-constant bulk-rock composition

In previous example, we have seen that excess water can trigger large amount of partial melting. One way to make the model slightly more realistic is by dynamically adjusting the composition at every step by removing excess water from it.

This can be be done by modifying your calculation as follow:

```julia
if "H2O" in out[i].ph
    id_h2o      = findfirst(out[i].ph .== "H2O")
    h2o_wt      = out[i].ph_frac_wt[id_h2o]
    h2o_comp_wt = out[i].PP_vec[id_h2o - out[i].n_SS].Comp_wt

    X           = X .- (h2o_wt .* h2o_comp_wt)
end
```
- Here, we first look for the id of "H2O" pure phase. Then, we get the water fraction in `wt`. Subsequently, we retrieve the composition of "H2O".  This part is slightly more complex as the information of solution models (`SS_vec`) are stored in a different substructure with respect to pure phases (`PP_vec`). In order to find the right id for "H2O" in the the `PP_vec`substructure we do `id_h2o - out[i].n_SS` where `out[i].n_SS` is the total number of solution models in the stable assemblage.

> [!NOTE]
> The previous code snipped has to be placed after calling `single_point_minimization()`

> Mind that for the igneous database, there is a fluid model "fl" instead of pure water.



<!-- 
```julia
data        =   Initialize_MAGEMin("mp", verbose=-1);
geo = 10.0
z   = collect(range(10.0,100,n))
P   = z .* 1000.0 .* 3300 .* 9.81 ./ 1e8
T   = 100.0 .+ z .* geo
Xup = copy(X)
for i=1:n
    T_calc  = T[i]       # retrieves the temperature from the temperature array we just defined  
    P_calc  = P[i]
    out[i]  = single_point_minimization(P_calc, T_calc, data, X=Xup, Xoxides=Xoxides, sys_in=sys_unit)

    if "H2O" in out[i].ph
        id_h2o      = findfirst(out[i].ph .== "H2O")
        h2o_wt      = out[i].ph_frac_wt[id_h2o]
        h2o_comp_wt = out[i].PP_vec[id_h2o - out[i].n_SS].Comp_wt

        Xup         = Xup .- (h2o_wt .* h2o_comp_wt)
    end
end

```
-->

