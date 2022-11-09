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

Application
    ausmPlusUpFlux

Description
    AUSM+up flux class, computes numerical fluxes based on AUSM+ scheme. Uses primitive 
    variables as input and gives back conservative numerical fluxes.

    Meng-Sing Liou, "A sequel to AUSM, PartII: AUSM+ -up for all speeds"
    Journal of Computational Physics 214 (2006): 137-170

Author
    Sebastian Saegeler  All rights reserved.
    Rewritten by Jiri Furst

SourceFiles
    AUSMplusUpFlux.H

\*---------------------------------------------------------------------------*/

#ifndef AUSMplusUpFlux_H
#define AUSMplusUpFlux_H

#include "dbnsFlux.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class AUSMplusUpFlux Declaration
\*---------------------------------------------------------------------------*/

class AUSMplusUpFlux: public dbnsFlux
{
private:
    AUSMplusUpFlux();
    AUSMplusUpFlux(const AUSMplusUpFlux&);
  
public:
    TypeName("AUSMplusUpFlux");

    // Constructors

    //- Constructor 
    AUSMplusUpFlux(const fvMesh&, const fluidThermo& thermo, const dictionary&);


    // Destructor - default
    virtual ~AUSMplusUpFlux() {};

    // Member Functions
    
    //- Update fluxes based on current left and right states
    virtual void evaluateFlux
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
    ) const;

private:
  
    scalar beta_;
    scalar MaInf_;
    scalar Kp_;
    scalar Ku_;
    scalar sigma_;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
