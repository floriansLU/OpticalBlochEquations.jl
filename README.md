# OpticalBlochEquations

[![Build Status](https://github.com/floriansLU/OpticalBlochEquations.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/floriansLU/OpticalBlochEquations.jl/actions/workflows/CI.yml?query=branch%3Amaster)

## Basic functionality
The package OpticalBlochEquations.jl is suitable for use in a REPL or in a Jupyter notebook. 

A typical session for producing the output similar to Fig.~\ref{fig:doppler} is shown in the following example. 

First, the necessary packages must be loaded. 

```julia
using OpticalBlochEquations
using Test
using Plots
using DelimitedFiles
gr()
```

Then, basic parameters of the transition and laser radiation must be defined.

```julia
params = param(cesiumD1)
laser_params = laser()
```

Then, it is necessary to define the polarization vectors for excitation, observation, and probing:

```julia
e_vec_ex = ElectricVector(1, π / 2, 0).cyclic
e_vec_obs = ElectricVector(1, π / 2, π / 2).cyclic
e_vec_probe = ElectricVector(0, π / 2, π / 4).cyclic
evecs=(e_vec_ex,e_vec_obs,e_vec_probe)
```

Finally, we define the number of steps in the average over the Doppler profile and the magnetic field:

```julia
Doppler_steps=150
B₀=0
```

Now we can compute the signals. The output of the function is a tuple that contains the fluorescence, the absorption, the ground-state density matrix, and the excited-state density matrix.

```julia
signals(B₀, params, laser_params, evecs, Doppler_steps)
```

## Plotting a probability surface
We can plot a probability surface of the atomic angular momentum distribution under the calculated conditions (magnetic field, laser, polarization) as follows
```julia
    res=signals_for_pmap(B₀)
    ρgg=res[3]
    ρee=res[4]
    ρggJ4=ρgg[1:9,1:9]
    J=4
    surface=plotProbSurf(J,ρggJ4)
    show(surface)
```


## Distributed Computing
In the event that we would like to calculate the signals over a wide range of magnetic-field values, it might be advantageous to use the Distributed.jl package with the function pmap. 
In that case, before loading the packages above, we should load the Distributed.jl package and start the macro @everywhere:

```julia
using Distributed
@everywhere begin
```

Then, after having defined all the parameters, we define a function for use with pmap:

```julia
signals_for_pmap(B₀) = signals(B₀, params, laser_params, evecs, Doppler_steps)
```

Now we can define a range of magnetic field values and run the program in multiple processes using pmap:

```julia
Brange=-40:5:40  
res = map(signals_for_pmap, Brange)
```

To extract the results, plot them and save them to a file, the following procedure can be used:

```julia
CsI=[]
CsA=[]

for i in 1:length(Brange)
    push!(CsA,res[1][i][2] |> real)
    push!(CsI,res[1][i][1] |> real)
end

default(legend=false)
p = plot(Brange, CsA)
plot!(p, title="A ar Doplera efektu", xlabel="B")
savefig(p, "A-ar-doplera-efektu.png")

p = plot(Brange, CsI)
plot!(p, title="I ar Doplera efektu", xlabel="B")
savefig(p, "I-ar-doplera-efektu.png")

open("I_Doplera.txt", "w") do io       
    writedlm(io, [Brange CsI] )               
end                                   

open("A_Doplera.txt", "w") do io       
    writedlm(io, [Brange CsA] )                
end
```
'''
