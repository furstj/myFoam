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

#include "AUSMPWplusFlux.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(AUSMPWplusFlux, 0);
    addToRunTimeSelectionTable(dbnsFlux, AUSMPWplusFlux, dictionary);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::AUSMPWplusFlux::AUSMPWplusFlux(const fvMesh&, const fluidThermo& thermo, const dictionary& dict):
    dbnsFlux(thermo)
{
    dictionary mySubDict( dict.subOrEmptyDict("AUSMPWplusFluxCoeffs") );
    alpha_ = mySubDict.lookupOrAddDefault("alpha", 3.0/16.0);
    beta_  = mySubDict.lookupOrAddDefault("beta", 1.0/8.0);
    
    if (mySubDict.lookupOrDefault("printCoeffs", false))
        Info << mySubDict << nl;
};


void Foam::AUSMPWplusFlux::evaluateFlux
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
    const scalar& RLeft,
    const scalar& RRight,
    const scalar& CvLeft,
    const scalar& CvRight,
    const vector& Sf,
    const scalar& magSf,
    const scalar& meshPhi
) const
{
    // Step 1: decode left and right:
    // normal vector
    const vector normalVector = Sf/magSf;
    
    const scalar qMesh = meshPhi / magSf;

    // Ratio of specific heat capacities
    const scalar kappaLeft = (RLeft + CvLeft)/CvLeft;
    const scalar kappaRight = (RRight + CvRight)/CvRight;
    const scalar kappa = 0.5 * (kappaLeft + kappaRight);

    // Density
    const scalar rhoLeft = pLeft/(RLeft*TLeft);
    const scalar rhoRight = pRight/(RRight*TRight);

    // Normal velocity (including mesh movement)
    const scalar qLeft = (ULeft & normalVector) - qMesh;
    const scalar qRight = (URight & normalVector) - qMesh;
    
    // "normal" total enthalpy
    const scalar HnLeft  = CvLeft*TLeft + pLeft/rhoLeft + 0.5*sqr(qLeft);
    const scalar HnRight = CvRight*TRight + pRight/rhoRight + 0.5*sqr(qRight);
    const scalar Hn = (HnLeft + HnRight) / 2.0;
    
    const scalar aStar = sqrt(2.0*(kappa-1)/(kappa+1) * Hn);
    const scalar af = (qLeft +  qRight > 0 ?
                       sqr(aStar)/max(fabs(qLeft), aStar) :
                       sqr(aStar)/max(fabs(qRight), aStar) );

    const scalar MLeft  = qLeft / af;
    const scalar MRight = qRight / af;

    scalar Mlp, plp;
    if (MLeft >= 1.)
    {
        Mlp = MLeft;
        plp = 1.0;
    }
    else if (MLeft > -1.0)
    {
        Mlp = sqr(1 + MLeft) / 4 + beta_ * sqr(1 - sqr(MLeft));
        plp = sqr(1 + MLeft)*(2 - MLeft)/4 + alpha_*MLeft*sqr(1 - sqr(MLeft));
    }
    else
    {
        Mlp = 0.0;
        plp = 0.0;
    }
    
    scalar Mrm, prm;
    if (MRight <= -1.)
    {
        Mrm = MRight;
        prm = 1.0;
    }
    else if (MRight < 1.0)
    {
        Mrm = - sqr(1 - MRight) / 4 - beta_ * sqr(1 - sqr(MRight));
        prm = sqr(1 - MRight)*(2 + MRight)/4 - alpha_*MRight*sqr(1 - sqr(MRight));
    }
    else
    {
        Mrm = 0.0;
        prm = 0.0;
    }

    const scalar Mf = Mlp + Mrm;
    const scalar pf = plp*pLeft + prm*pRight;
    
    const scalar aw  = 1 - pow(min(pLeft/pRight, pRight/pLeft),3);

    const scalar afl = (pf > 0 ? pLeft/pf - 1 : 0.0);
    const scalar afr = (pf > 0 ? pRight/pf - 1 : 0.0);

    if (Mf > 0)
    {
        Mlp += Mrm*( 1 - aw*(1+afr) + (afr-afl) );
        Mrm *= aw * (1+afr);
    }
    else
    {
        Mrm += Mlp*( 1 - aw*(1+afl) + (afl-afr) );
        Mlp *= aw * (1+afl);
    }

    const scalar rhoHLeft  = rhoLeft*(CvLeft*TLeft + 0.5*magSqr(ULeft)) + pLeft;
    const scalar rhoHRight = rhoRight*(CvRight*TRight + 0.5*magSqr(URight)) + pRight;
    
    rhoFlux  = magSf*af*(Mlp*rhoLeft + Mrm*rhoRight);
    rhoUFlux = magSf*af*(Mlp*rhoLeft*ULeft + Mrm*rhoRight*URight) + pf*Sf;
    rhoEFlux = magSf*af*(Mlp*rhoHLeft + Mrm*rhoHRight);
    
}


// ************************************************************************* //
