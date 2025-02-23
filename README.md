# OpticalBlochEquations
[![Build Status](https://github.com/floriansLU/OpticalBlochEquations.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/floriansLU/OpticalBlochEquations.jl/actions/workflows/CI.yml?query=branch%3Amaster)

## Introduction
The Optical Bloch Equations (OBEs) are useful for calculating the evolution of the density matrix of an atomic ensemble under the action of some Hamiltonian. A common situation concerns atoms with hyperfine structure that interact with an external magnetic field and laser radiation. We present a toolkit for solving the OBEs based on the \textit{QuantumOptics.jl} package in the Julia language. Using these tools makes the code much more readable than previous implementations in C/C++, but almost as fast and easier to parallelize. The toolkit includes functions for calculating the steady-state solution of density matrix of alkali metal atoms in the presence of an external magnetic field and exposed to a pump laser beam of arbitrary polarization and propagation direction. Based on this density matrix, the toolkit offer functions to determine the fluorescence intensity of arbitrary polarization and direction as well as the absorption of a weak probe beam, also of arbitrary polarization and propagation direction. The density matrix can be obtained and from the density matrix, angular momentum distributions can be calculated and plotted.

## Activating the project
Before doing anything, it is necessary to install `Julia`, which can be downloaded from [the Julia website](https://julialang.org/). This program has been tested with the latest version of Julia, currently v1.11.3. The package manager should be able to handle older or newer versions.  It is also highly recommended to install the package [IJulia](https://github.com/JuliaLang/IJulia.jl), which allows Jupyter notebooks to use the `Julia` kernel. Instructions are give in this [web page](https://julialang.github.io/IJulia.jl/stable/manual/running/). The notebooks in the `examples/` folder have been tested using the [Anaconda Navigator](https://www.anaconda.com/). Instructions for installing Jupyter notebooks for Julia in Anaconda can be found [here](https://www.engineered-mind.com/coding/install-julia-jupyter-notebook/). Jupyter notebooks can also be used through [VS Code](https://code.visualstudio.com/docs/languages/julia), as shown in this [tutorial](https://www.matecdev.com/posts/julia-introduction-vscode.html) .   

At the moment, the package has not been registered yet. Therefore, to use it, it must be activated in your environment. Create a project directory and perform the following steps
```
shell> git clone https://github.com/floriansLU/OpticalBlochEquations.jl.git
Cloning into 'OpticalBlochEquations.jl'...
...

(@v1.11) pkg> activate OpticalBlochEquations.jl
Activating project at `~/OpticalBlochEquations.jl`

(Example) pkg> instantiate
  No Changes to `~/OpticalBlochEquations.jl/Project.toml`
  No Changes to `~/OpticalBlochEquations.jl/Manifest.toml`
```


## Basic functionality
The package OpticalBlochEquations.jl is suitable for use in a REPL or in a Jupyter notebook. 
A typical session is shown in the following example. 

First, the necessary packages must be loaded. 

```julia
using OpticalBlochEquations
using Test
using Plots
using DelimitedFiles
gr()
```

Then, basic parameters of the transition and laser radiation must be defined.

## Defining the parameters
Now we set up the calculation. The function `param(cesiumD1)` defines the quantum numbers of the ground and excited states of the Cesium D1 line and create a structure that contain all the necessary operators. The function laser creates a structure for the Rabi frequency, natural line width, transit relaxation, laser frequency, laser linewidth, temperature, atomic mass, and other constants.  



```julia
params = param(cesiumD1)
laser_params = laser()
```

Then, it is necessary to define the polarization vectors for excitation, observation, and probing.
The interation of the atomic dipole with the electric field of the light is given by: 

$\hat{V}=-\hat{\mathbf{d}}\cdot\hat{\mathbf{E}}(t)$

The electric field is given by

$\hat{\mathbf{E}}(t)=\varepsilon(t)\mathbf{\epsilon} + \varepsilon^{*}(t)\mathbf{\epsilon^{*}}$

where $\epsilon$ is the polarization vector and 

$\varepsilon(t)=| \varepsilon_{\overline{\omega}}|e^{i\Phi(t)-i(\overline{\omega}-\mathbf{k}_{\overline{\omega}}\cdot \mathbf{v})t}$,

where $\varepsilon_{\overline{\omega}}$ is the amplitude,$\Phi(t)$ is a, possibly, time-dependent phase, $\overline{\omega}$ is the frequency of the light, $\mathbf{k}_{\overline{\omega}}$ is the wave vector of the light, and $\mathbf{v}$ is the velocity of the atom that interacts with the light.    

The polarization is defined using spherical polarization vectors. The spherical polarization vectors correspond to light that is left-circularly polarized ($\epsilon^{+1}$), linearly polarized along the $z$-direction ($\epsilon^{0}$), or right-circularly polarized ($\epsilon^{-1}$):. 

$\epsilon^{+1} = -\frac{1}{\sqrt{2}}\left( \epsilon_x -i \epsilon_y \right)$

$\epsilon^{0} = \epsilon_z$

$\epsilon^{-1} = \frac{1}{\sqrt{2}}\left( \epsilon_x +i \epsilon_y \right)  $

where $\epsilon_q$, $q \in \{x,y,z\}$ are the polarization vector's projection onto the Cartesian coordinate axes. 

Thus, \[1,0,0] would correspond to left-circularly polarized light with electric field vector rotating in the $xy$-plane, \[0,1,0] to light that is linearly polarized along the $z$-axis, and \[0,0,1] would correspond to right-circularly polarized light with the electric-field vector rotating in the $xy$-plane.  

We can rotate into an arbitrary coordinate system using the Euler angles $\alpha$, $\beta$, and $\gamma$. Thus, for example, we could start with linearly polarized light \[0,1,0] and rotate it by $\beta=\pi /2$ into the $xy$-plane to obtain counter-rotating left- and right-circularly polarized light. 

For more details on the polarization, see the following:

\[1]  M. Auzinsh, D. Budker, S. Rochester, Optically Polarized Atoms: Understanding Light-atom Interactions, OUP Oxford, 2010.

\[2] D. A. Varshalovich, A. N. Moskalev, V. K. Khersonskii, [Quantum Theory of Angular Momentum](https://library.oapen.org/handle/20.500.12657/50493), World Scientific Co. Pte. Ltd., Singapore, 2011.


```julia
e_vec_ex = ElectricVector([1,0,0], 0, π / 2, 0).cyclic       # left-circularly polarized light; the rotation rotates the propagation from the $z$-direction to the $x$-direction.
e_vec_obs = ElectricVector([1,0,0], 0, π / 2, π / 2).cyclic  # left-circularly polarized light; the rotation rotates the propagation from the $z$-direction to the $x$-direction.
e_vec_probe = ElectricVector([0,1,0], 0, π / 2, π / 4).cyclic  # linearly polarized light light whose polarization vector is rotated from being parallel to the $z$-axis to being half-way between the $x$- and $y$-axis in the $xy$-plan.
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

More details can be found in the Jupyter notebooks ```examples/testing-OBE-D1-transitions.ipynb``` and ```examples/testing-OBE-Rabi-frequency.ipynb```

## Plotting a probability surface
 The probability of a certain angular momentum value can be expressed as a surface in 3-dimensional space as follows:

    $\rho_{FF}(\theta,\phi)=\sqrt{\frac{4\pi}{2J+1}}\sum_{\kappa=0}^{2F}\sum_{q=-\kappa}^{\kappa} \left< F F \kappa 0 | FF \right> \rho^{\kappa q}Y_{\kappa q}(\theta,\phi)$.


The polarization moments $\rho^{\kappa q}$ can be expressed as

```math
\rho^{\kappa q} =\mathrm{Tr}\left(\rho \mathcal{T}_{q}^{\kappa}\right) \\ 

=\sum_{mm'}\rho_{mm'}\left(\mathcal{T}_q^{\kappa} \right)_{mm'}  \\

=\sum_{mm'}(-1)^{F-m}\left< F m' F, -m| \kappa q\right>\rho_{mm'} ,
```

where the $` \mathcal{T}_q^{\kappa} `$ are irreducible tensor operators called polarization operators and the density matrix is defined as

```math
    \rho=\sum_{mm'}\rho_{mm'}\left| m \right> \left< m' \right|.
```

We can plot a probability surface of the atomic angular momentum distribution under the calculated conditions (magnetic field, laser, polarization) as follows
```julia
    res=signals_for_pmap(B₀)
    ρgg=res[3]
    ρee=res[4]
    ρggJ4=ρgg[1:9,1:9]
    J=4
    surface=plotProbSurf(J,ρggJ4)
```
More details can be found in the Jupyter notebook ```examples/Angular-Momentum-Surface.ipynb```

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
More details can be found in the Jupyter notebook ```examples/parallel-OBE.ipynb```

## Program output and plotting
The program returns a tuple containing the fluorescence of the pump beam of the specified polarization, the absorption of a weak probe beam of specified polarization, the density matrix of the ground state $\rho_{g_ig_j}$ and the density matrix of the excited state $\rho_{e_ie_j}$.  We can run the program over a a range of magnetic field values defined in `Brange`.


```julia
Brange=-40:5:40 
res=[]
for B₀ in Brange
    result=signals(B₀, par, laz, evecs, Doppler_steps)
    append!(res,result)
end
```

To extract the probe absorption, plot it as a function of magnetic field, and save the results to a file, the following procedure can be used:

```julia
x=Brange
myResults=reshape(res,(4,length(res)))
y=myResults[2,:] |> real;
plot(x,y,xlabel="Magnetic Field [G]",ylabel="Probe absorption")
```
More details can be found in the Jupyter notebooks ```examples/testing-OBE-D1-transitions.ipynb``` and ```examples/testing-OBE-Rabi-frequency.ipynb```


```
'''
