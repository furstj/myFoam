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

#include "rotatedFlux.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(rotatedFlux, 0);
    addToRunTimeSelectionTable(dbnsFlux, rotatedFlux, dictionary);
}

// * * * * * * * * * * * * * * * Constructor  * * * * * * * * * * * * * //
Foam::rotatedFlux::rotatedFlux(const fvMesh& mesh, const fluidThermo& thermo, const dictionary& dict):
    dbnsFlux(thermo)
{
    dictionary subDict( dict.subDict("rotatedFluxCoeffs") );
    epsilon_ = subDict.lookupOrAddDefault("epsilon", 1e-3);

    diffusiveFlux_ = Foam::dbnsFlux::New(mesh, thermo, subDict, "diffusiveFlux");
    lowDiffusionFlux_ = Foam::dbnsFlux::New(mesh, thermo, subDict, "lowDiffusionFlux");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rotatedFlux::evaluateFlux
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
  if (mag(meshPhi)>0.0) 
    {
      FatalError
        << "This dbnsFlux is not ready to run with moving meshes." << nl
        << exit(FatalError);
    };

  const vector deltaQ = URight - ULeft;
  const scalar magDeltaQ = mag(deltaQ);

  if (magDeltaQ > epsilon_)
  {
      scalar rhoFlux1, rhoFlux2, rhoEFlux1, rhoEFlux2;
      vector rhoUFlux1, rhoUFlux2;

      const vector n = Sf / magSf;
      
      vector n1 = deltaQ/magDeltaQ;
      if ( (n1 & n) < 0) n1 *= -1;
      
      diffusiveFlux_->evaluateFlux(rhoFlux1, rhoUFlux1, rhoEFlux1, 
      pLeft, pRight, ULeft, URight, TLeft, TRight, RLeft, RRight, CvLeft, CvRight,
      n1*magSf, magSf, meshPhi);

      vector n2 = (n1^n)^n1;
      n2 /= mag(n2);

      lowDiffusionFlux_->evaluateFlux(rhoFlux2, rhoUFlux2, rhoEFlux2, 
      pLeft, pRight, ULeft, URight, TLeft, TRight, RLeft, RRight, CvLeft, CvRight,
      n2*magSf, magSf, meshPhi);

      const scalar alpha1 = n & n1;
      const scalar alpha2 = n & n2;

      rhoFlux = alpha1*rhoFlux1 + alpha2*rhoFlux2;
      rhoUFlux = alpha1*rhoUFlux1 + alpha2*rhoUFlux2;
      rhoEFlux = alpha1*rhoEFlux1 + alpha2*rhoEFlux2;
  }
  else
  {
      // According to Nishikawa, n1 is othogonal to n (i.e. alpha1=0)
      lowDiffusionFlux_->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux, 
      pLeft, pRight, ULeft, URight, TLeft, TRight, RLeft, RRight, CvLeft, CvRight,
      Sf, magSf, meshPhi);
  }
  
}

// ************************************************************************* //
