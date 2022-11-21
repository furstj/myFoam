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

#include "totalPressureRG.H"
#include "fluidThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(totalPressureRG, 0);
    addToRunTimeSelectionTable(functionObject, totalPressureRG, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::functionObjects::totalPressureRG::pTot(
    scalar p, scalar T, vector U, Foam::gasProperties& gasProps
)
{
    scalar S = gasProps.S(p, T);
    scalar H = gasProps.Hs(p, T) + 0.5*magSqr(U);
    scalar T0 = T;
    scalar p0 = p;

    const scalar tol = 1.e-8;
    const label maxIter = 100;
    
    label iter = 0;
    scalar dp, dT;
    do
    {
        scalar Cp = gasProps.Cp(p0, T0);
        scalar beta_p = gasProps.beta_p(p0, T0);
        scalar v = 1.0/gasProps.rho(p0, T0);
        
        scalar dH = H - gasProps.Hs(p0, T0);
        scalar dS = S - gasProps.S(p0, T0);
        dT = dH/Cp;
        dp = (Cp/T0*dT - dS)/(v*beta_p);
        
        
        T0 += dT;
        p0 += dp;
        
        if (iter++ > maxIter)
        {
            FatalErrorInFunction
                << "Maximum number of iterations exceeded: " << maxIter
                    << " T  : " << T0
                    << " p  : " << p0
                    << " Z  : " << gasProps.Z(p0,T0)
                    << " Cp : " << gasProps.Cp(p0,T0)
                    << " tol: " << tol
                    << abort(FatalError);
        }
        
    } while ( (mag(dp) > p*tol) || (mag(dT) > T*tol) );
    
    return p0;
}


bool Foam::functionObjects::totalPressureRG::calc()
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
                dimensionedScalar("total(p)", dimPressure, 0.0)
            );

        volScalarField& result = tresult.ref();
        
        forAll(mesh_.cells(), i)
        {
            result[i] = this->pTot(p[i], T[i], U[i], gasProps);
        };

        forAll(result.boundaryField(), patchi)
        {
            auto& presult  = result.boundaryFieldRef()[patchi];
            const auto& pU = U.boundaryField()[patchi];
            const auto& pp = p.boundaryField()[patchi];
            const auto& pT = T.boundaryField()[patchi];
            forAll(presult, i)
            {
                presult[i] = this->pTot(pp[i], pT[i], pU[i], gasProps);
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

Foam::functionObjects::totalPressureRG::totalPressureRG
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "U")
{
    setResultName("total(p)", "U");
}


// ************************************************************************* //
