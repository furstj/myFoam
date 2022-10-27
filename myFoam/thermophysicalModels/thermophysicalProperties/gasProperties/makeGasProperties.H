/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::makeGasProperties

Description

SourceFiles
    makeGases.C

\*---------------------------------------------------------------------------*/

#ifndef makeGasProperties_H
#define makeGasProperties_H

#include "genericGasProperties.H"
#include "makeThermo.H"

namespace Foam
{

#define defineGasProperties(ThermoPhys)                                 \
                                                                        \
    typedef genericGasProperties<ThermoPhys> gas##ThermoPhys;           \
                                                                        \
    defineTemplateTypeNameAndDebugWithName                              \
    (                                                                   \
        gas##ThermoPhys,                                                \
        ("gasProperties<" + ThermoPhys::typeName() + ">").c_str(),      \
        0                                                               \
    ); 

#define addGasProperties(Properties)                                    \
                                                                        \
    addToRunTimeSelectionTable                                          \
    (                                                                   \
    gasProperties,                                                      \
    Properties,                                                         \
    dictionary                                                          \
    );


#define makeGasProperties(Transport,Type,Thermo,EqnOfState,Specie)      \
                                                                        \
    typedefThermoPhysics(Transport,Type,Thermo,EqnOfState,Specie);      \
                                                                        \
    defineGasProperties(Transport##Type##Thermo##EqnOfState##Specie);   \
                                                                        \
    addGasProperties(gas##Transport##Type##Thermo##EqnOfState##Specie);
}    

#endif

// ************************************************************************* //
