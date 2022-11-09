/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "hllLMFlux.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(hllLMFlux, 0);
    addToRunTimeSelectionTable(dbnsFlux, hllLMFlux, dictionary);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hllLMFlux::evaluateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const scalar& pLeft,
    const scalar& pRight,
    const vector& ULeft,
    const vector& URight,
    const scalar& TLeft,
    const scalar& TRight,
    const vector& Sf,
    const scalar& magSf,
    const scalar& meshPhi
) const
{
    //if (mag(meshPhi)>0.0) 
    //{
    //  FatalError
    //    << "This dbnsFlux is not ready to run with moving meshes." << nl
    //    << exit(FatalError);
    //};

    // Step 1: decode left and right:
    // normal vector
    const vector normalVector = Sf/magSf;
    
    const scalar qMesh = meshPhi / magSf;

    // Compute conservative variables 

    // Density
    const scalar rhoLeft = gas().rho(pLeft, TLeft);
    const scalar rhoRight = gas().rho(pRight, TRight);

    // DensityVelocity
    const vector rhoULeft = rhoLeft*ULeft;
    const vector rhoURight = rhoRight*URight;

    // Compute left and right total enthalpies:
    const scalar HLeft = gas().Hs(pLeft, TLeft) + 0.5*magSqr(ULeft);
    const scalar HRight = gas().Hs(pRight, TRight) + 0.5*magSqr(URight);

    // DensityTotalEnergy
    const scalar rhoELeft = rhoLeft*HLeft - pLeft;
    const scalar rhoERight = rhoRight*HRight - pRight;


    // Compute qLeft and qRight (q_{l,r} = U_{l,r} \bullet n)
    const scalar qLeft = (ULeft & normalVector) - qMesh;
    const scalar qRight = (URight & normalVector) - qMesh;

    // Speed of sound, for left and right side
    const scalar aLeft  = Foam::sqrt(max(0.0, gas().c(pLeft, TLeft)));
    const scalar aRight = Foam::sqrt(max(0.0, gas().c(pRight, TRight)));

    // Step 2:
    // needs rho_{l,r}, U_{l,r}, H_{l,r}, kappa_{l,r}, Gamma_{l,r}, q_{l,r}

    // Compute Roe weights
    const scalar rhoLeftSqrt = Foam::sqrt(max(0.0,rhoLeft));
    const scalar rhoRightSqrt = Foam::sqrt(max(0.0,rhoRight));

    const scalar wLeft = rhoLeftSqrt
        /stabilise((rhoLeftSqrt + rhoRightSqrt),VSMALL);

    const scalar wRight = 1 - wLeft;

    // Roe averaged velocity
    const vector UTilde = wLeft*ULeft + wRight*URight;

    // Roe averaged contravariant velocity
    const scalar contrUTilde = (UTilde & normalVector);

    // Roe averaged total enthalpy
    const scalar HTilde = wLeft*HLeft + wRight*HRight;

    // Speed of sound with Roe reconstruction values
    const scalar pTilde = wLeft*pLeft + wRight*pRight;
    const scalar TTilde = gas().THs(HTilde - 0.5*sqr(contrUTilde), pTilde, (TLeft+TRight)/2);
    const scalar aTilde = gas().c(pTilde, TTilde);

    // Step 3: compute signal speeds for face:
    const scalar SLeft  = min(qLeft - aLeft, contrUTilde - qMesh - aTilde);
    const scalar SRight = max(contrUTilde - qMesh + aTilde, qRight + aRight);



    if (pos(SLeft))
    {
        rhoFlux  = (qLeft*rhoLeft)*magSf;
        rhoUFlux = (qLeft*rhoULeft + pLeft*normalVector)*magSf;
        rhoEFlux = (qLeft*(rhoELeft + pLeft) + pLeft*qMesh)*magSf;
    }
    else if (neg(SRight))
    {
        rhoFlux  = (qRight*rhoRight)*magSf;
        rhoUFlux = (qRight*rhoURight + pRight*normalVector)*magSf;
        rhoEFlux = (qRight*(rhoERight + pRight) + pRight*qMesh)*magSf;
    }
    else
    {
        rhoFlux  = (SRight*(qLeft*rhoLeft) - SLeft*(qRight*rhoRight))*magSf/(SRight - SLeft);
        rhoUFlux = (SRight*(qLeft*rhoULeft + pLeft*normalVector) -
                    SLeft*(qRight*rhoURight + pRight*normalVector))*magSf/(SRight - SLeft);
        rhoEFlux = (SRight*(qLeft*(rhoELeft + pLeft) + pLeft*qMesh) -
                    SLeft*(qRight*(rhoERight + pRight) + pRight*qMesh))*magSf/(SRight - SLeft);

        // Compute Mach number function and adjusted pressure in star region
        const scalar MachLeft = mag(ULeft)/aLeft;
        const scalar MachRight = mag(URight)/aRight;
        const scalar LIM = min(max(MachLeft,MachRight),1.0);

        rhoFlux += SRight*SLeft/(SRight - SLeft)*(rhoRight - rhoLeft)*magSf;
        rhoUFlux += LIM*SRight*SLeft/(SRight - SLeft)*(rhoURight - rhoULeft)*magSf;
        rhoEFlux += SRight*SLeft/(SRight - SLeft)*(rhoERight - rhoELeft)*magSf;
    }


}


// ************************************************************************* //
