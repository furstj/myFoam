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
    hllcLMFlux

Description
    HLLC flux class, computes numerical fluxes based on HLLC scheme. Uses
    primitive variables as input and gives back conservative numerical fluxes.

    Toro, E.I. and Spruce, M. and Spears, W.: Restoration of the Contact
    Surface in the HLL-Riemann Solver, Shock Waves, 4:25-34, 1996

	 Modified to give correct scaling of pressure and density fluctuations in
	 the limit of zero Mach number.

	 Xie, W. and Zhang, R. and Lai, J. and Li, H.:An accurate and robust 
	 HLLC-type Riemann solver for the compressible Euler system at various 
	 Mach numbers, Int. J. Numer. Meth. Fluids, 1-34, 2018

Author
    Oliver Borm
    Aleksandar Jemcov
    Rewrite by Hrvoje Jasak
	 Low Mach correction by Michael Groom

SourceFiles
    hllcLMFlux.H

\*---------------------------------------------------------------------------*/

#ifndef hllcLMFlux_H
#define hllcLMFlux_H

#include "dbnsFlux.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class hllcLMFlux Declaration
\*---------------------------------------------------------------------------*/

class hllcLMFlux: public dbnsFlux
{
public:
    TypeName("hllcLMFlux");

    // Constructors

    //- Constructor null
    hllcLMFlux(const fvMesh&, const fluidThermo& thermo, const dictionary&):
        dbnsFlux(thermo)
        {};


    // Destructor - default
    virtual ~hllcLMFlux() {};

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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
