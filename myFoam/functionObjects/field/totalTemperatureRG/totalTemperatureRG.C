/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "totalTemperatureRG.H"
#include "fluidThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(totalTemperatureRG, 0);
    addToRunTimeSelectionTable(functionObject, totalTemperatureRG, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::totalTemperatureRG::calc()
{
    if
    (
        foundObject<volVectorField>(fieldName_)
     && foundObject<fluidThermo>(fluidThermo::dictName)
    )
    {
        const fluidThermo& thermo =
            lookupObject<fluidThermo>(fluidThermo::dictName);

        autoPtr<gasProperties> pGasProps(gasProperties::New(thermo));
        gasProperties& gasProps = pGasProps.ref();
        
        const volVectorField& U = lookupObject<volVectorField>(fieldName_);
        const volScalarField& p = thermo.p();
        const volScalarField& T = thermo.T();
        const volScalarField  H = thermo.he() + 0.5*magSqr(U);

        auto tresult =
            tmp<volScalarField>::New
            (
                IOobject
                (
                    resultName_,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ
                ),
                mesh_,
                dimensionedScalar("total(T)", dimTemperature, 0.0)
            );

        volScalarField& result = tresult.ref();

        forAll(mesh_.cells(), i)
        {
            scalar s    = gasProps.S(p[i], T[i]);
            scalar pTot = gasProps.pHS(H[i], s, p[i]);
            result[i] = gasProps.THs(H[i], pTot, T[i]);
        };

        forAll(result.boundaryField(), patchi)
        {
            auto& presult  = result.boundaryFieldRef()[patchi];
            const auto& pH = H.boundaryField()[patchi];
            const auto& pp = p.boundaryField()[patchi];
            const auto& pT = T.boundaryField()[patchi];
            forAll(presult, i)
            {
                scalar s    = gasProps.S(pp[i], pT[i]);
                scalar pTot = gasProps.pHS(pH[i], s, pp[i]);
                presult[i] = gasProps.THs(pH[i], pTot, pT[i]);
            }
        }

        
        return store(
            resultName_,
            tresult
        );
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::totalTemperatureRG::totalTemperatureRG
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "U")
{
    setResultName("total(T)", "U");
}


// ************************************************************************* //
