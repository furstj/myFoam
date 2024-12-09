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

#include "AUSMplusUpFlux.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(AUSMplusUpFlux, 0);
    addToRunTimeSelectionTable(dbnsFlux, AUSMplusUpFlux, dictionary);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::AUSMplusUpFlux::AUSMplusUpFlux(const fvMesh&, const fluidThermo& thermo, const dictionary& dict):
    dbnsFlux(thermo)
{
    dictionary mySubDict( dict.subOrEmptyDict("AUSMplusUpFluxCoeffs") );
    beta_    = mySubDict.lookupOrAddDefault("beta", 1.0/8.0);
    MaInf_   = mySubDict.lookupOrAddDefault("MaInf", 0.3);
    Kp_      = mySubDict.lookupOrAddDefault("Kp", 0.25);
    Ku_      = mySubDict.lookupOrAddDefault("Ku", 0.75);
    sigma_   = mySubDict.lookupOrAddDefault("sigma", 1.0);
    deltaEF_ = mySubDict.lookupOrAddDefault("deltaEF", 1e-4);
    
    if (mySubDict.lookupOrDefault("printCoeffs", false))
        Info << mySubDict << nl;
};


void Foam::AUSMplusUpFlux::evaluateFlux
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
    // Step 1: decode left and right:
    // normal vector
    const vector normalVector = Sf/magSf;
    
    const scalar qMesh = meshPhi / magSf;

    // Density
    const scalar rhoLeft = gas().rho(pLeft, TLeft);
    const scalar rhoRight = gas().rho(pRight, TRight);

    // DensityTotalEnthalpy
    const scalar rhoHLeft  = rhoLeft*(gas().Hs(pLeft, TLeft) + 0.5*magSqr(ULeft));
    const scalar rhoHRight = rhoRight*(gas().Hs(pRight, TRight) + 0.5*magSqr(URight));

    // DensityVelocity
    const vector rhoULeft = rhoLeft*ULeft;
    const vector rhoURight = rhoRight*URight;

    // Compute qLeft and qRight (q_{l,r} = U_{l,r} \bullet n)
    const scalar qLeft = (ULeft & normalVector) - qMesh;
    const scalar qRight = (URight & normalVector) - qMesh;

    const scalar Ht = 0.5*(rhoHLeft/rhoLeft + rhoHRight/rhoRight);

    // Compute sonic parameters by solving h(p,T) + 0.5c^2(p,T) = H
    const scalar pStar = (pLeft + pRight)/2;
    scalar TStar = (TLeft + TRight)/2;
    for (label iter=0; iter<2; iter++)  // TODO: more iterations for real gas!
    {
        scalar Cp = gas().Cp(pStar, TStar);
        scalar Cv = Cp - gas().CpMCv(pStar, TStar);
        scalar f = gas().Hs(pStar, TStar) + 0.5*sqr(gas().c(pStar,TStar)) - Ht;
        scalar fPrime = 0.5*(Cp/Cv+1)*Cp;
        TStar -= f/fPrime;
    }
    
    const scalar aStar = gas().c(pStar, TStar);
    const scalar aHatLeft  = sqr(aStar) / max(aStar, qLeft);
    const scalar aHatRight = sqr(aStar) / max(aStar, -qRight);
    const scalar aTilde = min(aHatLeft, aHatRight);
    const scalar rhoTilde = 0.5*(rhoLeft+rhoRight);

    const scalar sqrMaDash = ((sqr(qLeft)+sqr(qRight))/(2.0*sqr(aTilde))+deltaEF_)/(1+deltaEF_);
    const scalar sqrMaZero = min(1.0,max(sqrMaDash,sqr(MaInf_)));
    const scalar MaZero    = Foam::sqrt(sqrMaZero);

    const scalar fa = MaZero*(2.0-MaZero);

    const scalar alpha = 3.0/16.0*(-4.0+5.0*sqr(fa));
    
    const scalar MaRelLeft  = qLeft /aTilde;
    const scalar MaRelRight = qRight/aTilde;

    const scalar magMaRelLeft  = (mag(MaRelLeft) +deltaEF_)/(1+deltaEF_);
    const scalar magMaRelRight = (mag(MaRelRight)+deltaEF_)/(1+deltaEF_);
    
    const scalar Ma1PlusLeft   = 0.5*(MaRelLeft +magMaRelLeft );
    const scalar Ma1MinusRight = 0.5*(MaRelRight-magMaRelRight);
    
    const scalar Ma2PlusLeft   =  0.25*sqr(MaRelLeft +1.0);
    const scalar Ma2PlusRight  =  0.25*sqr(MaRelRight+1.0);
    const scalar Ma2MinusLeft  = -0.25*sqr(MaRelLeft -1.0);
    const scalar Ma2MinusRight = -0.25*sqr(MaRelRight-1.0);
    
    const scalar Ma4BetaPlusLeft   = ((magMaRelLeft  >= 1.0) ? Ma1PlusLeft   : (Ma2PlusLeft  *(1.0-16.0*beta_*Ma2MinusLeft)));
    const scalar Ma4BetaMinusRight = ((magMaRelRight >= 1.0) ? Ma1MinusRight : (Ma2MinusRight*(1.0+16.0*beta_*Ma2PlusRight)));
        
    const scalar Mp = -Kp_/fa*max(1.0-sigma_*sqrMaDash,0.0)*(pRight-pLeft)/(rhoTilde*sqr(aTilde));

    const scalar P5alphaPlusLeft   = ((magMaRelLeft  >= 1.0) ?
    (Ma1PlusLeft/MaRelLeft)    : (Ma2PlusLeft  *(( 2.0-MaRelLeft) -16.0*alpha*MaRelLeft *Ma2MinusLeft )));
    const scalar P5alphaMinusRight = ((magMaRelRight >= 1.0) ?
    (Ma1MinusRight/MaRelRight) : (Ma2MinusRight*((-2.0-MaRelRight)+16.0*alpha*MaRelRight*Ma2PlusRight)));
    
    const scalar pU = -Ku_*P5alphaPlusLeft*P5alphaMinusRight*(rhoLeft+rhoRight)*(fa*aTilde)*(qRight-qLeft);
    
    const scalar MaRelTilde = Ma4BetaPlusLeft + Ma4BetaMinusRight + Mp;
    const scalar pTilde = pLeft*P5alphaPlusLeft + pRight*P5alphaMinusRight + pU;
    
    const scalar URelTilde = MaRelTilde*aTilde;
    const scalar magURelTilde = mag(MaRelTilde)*aTilde;
    // There is a typo in Luo et. al, J. Comp. Physics 194 (2004), Chap 4.2 Eq. 4.8
    // refer to the origial Paper from Liou, J. Comp. Physics 129 (1996), Chap4, Eq. 42
    rhoFlux  = (0.5*(URelTilde*(rhoLeft +rhoRight) -magURelTilde*(rhoRight -rhoLeft)))*magSf;
    rhoUFlux = (0.5*(URelTilde*(rhoULeft+rhoURight)-magURelTilde*(rhoURight-rhoULeft))+pTilde*normalVector)*magSf;
    rhoEFlux = (0.5*(URelTilde*(rhoHLeft+rhoHRight)-magURelTilde*(rhoHRight-rhoHLeft)) + pTilde*qMesh)*magSf;

}


// ************************************************************************* //
