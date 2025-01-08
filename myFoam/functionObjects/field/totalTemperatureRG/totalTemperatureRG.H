/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::functionObjects::totalEnthalpyRG

Group
    grpFieldFunctionObjects

Description
    Computes the total temperature as a \c volScalarField.

    Operands:
    \table
      Operand          | Type           | Location
      input            | volVectorField | $FOAM_CASE/\<time\>/\<inpField\>
      output file      | -              | -
      output field     | volScalarField | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    totalTemperatureRG1
    {
        // Mandatory entries (unmodifiable)
        type            totalTemperatureRG;
        libs            (fieldFunctionObjects);

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                          | Type | Req'd | Dflt
      type         | Type name: totalTemperatureRG           | word |  yes  | -
      libs         | Library name: myfieldFunctionObjects | word |  yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::fieldExpression
    - ExtendedCodeGuide::functionObjects::field::totalTemperatureRG

SourceFiles
    totalTemperatureRG.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_totalTemperatureRG_H
#define functionObjects_totalTemperatureRG_H

#include "fieldExpression.H"
#include "gasProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class totalTemperatureRG Declaration
\*---------------------------------------------------------------------------*/

class totalTemperatureRG
:
    public fieldExpression
{
    // Private Member Functions

    //- Calculate the total enthalpy field and return true if successful
    virtual bool calc();
    
public:

    //- Runtime type information
    TypeName("totalTemperatureRG");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        totalTemperatureRG
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        totalTemperatureRG(const totalTemperatureRG&) = delete;

        //- No copy assignment
        void operator=(const totalTemperatureRG&) = delete;


    //- Destructor
    virtual ~totalTemperatureRG() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
