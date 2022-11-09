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

#include "hllcLMFlux.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(hllcLMFlux, 0);
    addToRunTimeSelectionTable(dbnsFlux, hllcLMFlux, dictionary);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hllcLMFlux::evaluateFlux
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
  //  {
  //    FatalError
  //      << "This dbnsFlux is not ready to run with moving meshes." << nl
  //      << exit(FatalError);
  //  };

    // Step 1: decode left and right:
    // normal vector
    const vector normalVector = Sf/magSf;
    
    const scalar qMesh = meshPhi / magSf;

    // Compute conservative variables assuming perfect gas law

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

    // Speed of sound, for left and right side, assuming perfect gas
    const scalar aLeft = gas().c(pLeft, TLeft);    
    const scalar aRight = gas().c(pRight, TRight);


    const scalar MachLeft = mag(ULeft)/aLeft;
    const scalar MachRight = mag(URight)/aRight;

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

    // Static enthalpy using speed in normal direction
    const scalar hTilde = HTilde - 0.5*sqr(contrUTilde);

    // Roe averaged pressure (?)
    const scalar pTilde = wLeft*pLeft + wRight*pRight;

    const scalar TTilde = gas().THs(hTilde, pTilde, (TLeft+TRight)/2);
    
    const scalar aTilde = gas().c(pTilde, TTilde);

    // Step 3: compute signal speeds for face:
    const scalar SLeft  = min(qLeft - aLeft, contrUTilde - qMesh - aTilde);
    const scalar SRight = max(contrUTilde - qMesh + aTilde, qRight + aRight);

    const scalar SStar = (rhoRight*qRight*(SRight-qRight)
      - rhoLeft*qLeft*(SLeft - qLeft) + pLeft - pRight )/
        stabilise((rhoRight*(SRight-qRight)-rhoLeft*(SLeft-qLeft)),VSMALL);

    // Compute pressure in star region from the right side
    const scalar pStarRight =
        rhoRight*(qRight - SRight)*(qRight - SStar) + pRight;

    // Should be equal to the left side
    const scalar pStarLeft  =
        rhoLeft*(qLeft -  SLeft)*(qLeft - SStar) + pLeft;

    // Give a warning if this is not the case
    if (mag(pStarRight - pStarLeft) > 1e-6)
    {
        Info << "mag(pStarRight-pStarLeft) > VSMALL " << endl;
    }

    // Use pStarRight for pStar, as in theory, pStarRight == pStarLeft
    const scalar pStar = pStarRight;

	 // Compute Mach number function and adjusted pressure in star region
	 const scalar LIM = min(max(MachLeft,MachRight),1.0);
	 const scalar pStar2 = LIM*pStar+(1.0-LIM)*0.5*(pLeft+pRight);

    // Step 4: upwinding - compute states:
    scalar convectionSpeed = 0.0;
    scalar rhoState = 0.0;
    vector rhoUState = vector::zero;
    scalar rhoEState = 0.0;
    scalar pState = 0.0;
    scalar pState2 = 0.0;

    if (pos(SLeft))
    {
        // compute F_l
        convectionSpeed = qLeft;
        rhoState  = rhoLeft;
        rhoUState = rhoULeft;
        rhoEState = rhoELeft;
        pState = pLeft;
        pState2 = pLeft;
    }
    else if (pos(SStar))
    {
        scalar omegaLeft = scalar(1.0)/stabilise((SLeft - SStar), VSMALL);

        // Compute left star region
        convectionSpeed = SStar;
        rhoState  = omegaLeft*(SLeft - qLeft)*rhoLeft;
        rhoUState = omegaLeft*((SLeft - qLeft)*rhoULeft
        + (pStar2 - pLeft)*normalVector);
        rhoEState = omegaLeft*((SLeft - qLeft)*rhoELeft
        - pLeft*qLeft + pStar*SStar);
        pState = pStar;
        pState2 = pStar2;
    }
    else if (pos(SRight))
    {
        scalar omegaRight = scalar(1.0)/stabilise((SRight - SStar), VSMALL);

        // compute right star region
        convectionSpeed = SStar;
        rhoState  = omegaRight*(SRight - qRight)*rhoRight;
        rhoUState = omegaRight*((SRight - qRight)*rhoURight
        + (pStar2 - pRight)*normalVector);
        rhoEState = omegaRight*((SRight - qRight)*rhoERight
        - pRight*qRight + pStar*SStar);
        pState = pStar;
        pState2 = pStar2;
    }
    else if (neg(SRight))
    {
        // compute F_r
        convectionSpeed = qRight;
        rhoState  = rhoRight;
        rhoUState = rhoURight;
        rhoEState = rhoERight;
        pState = pRight;
        pState2 = pRight;
    }
    else
    {
        Info << "Error in HLLC Riemann solver" << endl;
    }

    rhoFlux  = (convectionSpeed*rhoState)*magSf;
    rhoUFlux = (convectionSpeed*rhoUState+pState2*normalVector)*magSf;
    rhoEFlux = (convectionSpeed*(rhoEState+pState) + pState*qMesh)*magSf;
}


// ************************************************************************* //
