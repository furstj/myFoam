# myFoam 

The repository contains several additional solvers for OpenFOAM. There are
several branches supporting different versions of OpenFOAM:

- **master** development branch supporting recent OpenFOAM (Foundation and ESI)
- **realGas** development branch with extended support of real gas
- **OF6** branch for OpenFOAM 6 and OpenFOAM v1812
- **OF5** branch for OpenFOAM 5


The package consists of two sub-directories:

- **myFoam** containing several extensions to standard OpenFOAM including:
  - some specific boundary conditions,
  - mesh movement solver,
  - Venkatakrisnan limiter.
- **Solver** containing
  - *myLusgsFoam*: the LU-SGS solver for steady or transient simulations
    of turbulent compressible flows, 
  - *mySonicLiquidFoam*: transient solver for trans-sonic/supersonic, laminar
    flow of a compressible liquid. 

Obsolete / removed parts:

  - *myRhoSimpleFoam*: steady state compressible gas with corrections for MRF, the OpenFOAM v2312 
     already has corrected version of EEqn. Therefore *myRhoSimpleFoam* was removed.


# Installation
The installation is quite easy, just do

    cd myFoam
    ./Allwmake
	cd ../mySolvers
	./Allwmake

## Real gas via CoolProp library
The **realGas** branch optionally supports also *CoolProp* library. In order to compile the package with*CoolProp* one has to set *COOLPROP* environment variable before compilation, i.e.

    cd myFoam
    export COOLPROP=/path/to/coolprop
    ./Allwmake

The compilation then looks for includes in *$COOLPROP/include* and for shared library in *$COOLPROP/lib*.

Note that the support for *CoolProp* is stiil highly experimental and the calculation is deadly slow when comparing with ideal gas or with cubic equation of state.


# Tutorial case
There is simple tutorial case in *tutorials/myLuSgsFoam*. Just run `./Allrun`
inside that directory.

# Notes
- **Compatibility with OpenFOAM v2012** there is a change in thermodynamics.
  The sensible enthalpy of ideal gas is calculated in v2012 as `h = cp*(T - Tref) + href`
  whereas `h = cp*T` in eralier versions. Therefore the solver requires to set
  `Tref = 0` in `thermophysicalProperties`, see tutorials.
  
# How to cite
Please cite as

    @article{Furst2018,
    author = {F{\"{u}}rst, Ji\v{r}{\'{\i}}},
    doi = {10.1016/j.compfluid.2018.04.020},
    issn = {00457930},
    journal = {Computers {\&} Fluids},
    month = {aug},
    pages = {332--339},
    title = {Development of a coupled matrix-free {LU-SGS} solver for turbulent compressible flows},
    volume = {172},
    year = {2018}
	}

	
