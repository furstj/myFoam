# myFoam 

The repository contains several additional solvers for OpenFOAM. There are
several branches supporting different versions of OpenFOAM:

- **master** development branch supporting recent OpenFOAM (Foundation and ESI)
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

# Installation
The installation is quite easy, just do

    cd myFoam
    ./Allwmake
	cd ../mySolvers
	./Allwmake

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

	
