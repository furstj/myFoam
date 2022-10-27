# myThermophysicalModels
Library of additional thermophysical models for OpenFOAM.

## Additional thermophysical models
* Aungier-Redlich-Kwong EOS

## gasProperties
The library implements a **gasProperties** library (modification of **thermophysicalProperties**) which allows the access to low-level thermophysical library (e.g. *rho(p,T)*). The gas properties are constructed from the same data as the standard thermo package using *constant/thermophysicalProperties*.

There is an example in *tutorial* folder.
