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

Class
    rotatedFlux


Description
    Rotated flux is a hybrid flux using a diffusive flux in shock-normal direction 
    a low diffusion scheme in othe directions. See

    H. Nishikawa, K. Kitamura, 
    “Very simple, carbuncle free, boundary-layer-resolving, rotated-hybrid Riemann solvers,” 
    Journal of Computational Physics, vol. 227, pp. 2560–2581, 2008.

    Corresponding part of fvSchemes should be:

    dbnsFlux rotatedFlux;

    rotatedFluxCoeffs
    {
        epsilon 1.e-3;
        diffusiveFlux rusanovFlux;
        lowDiffusionFlux roeFlux;
    }

    where epsilon serves as a stabilization parameter in calculation of n1

Author
    Jiri Furst

SourceFiles
    rotatedFlux.H

\*---------------------------------------------------------------------------*/

#ifndef rotatedFlux_H
#define rotatedFlux_H

#include "dbnsFlux.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class rotatedFlux Declaration
\*---------------------------------------------------------------------------*/

class rotatedFlux:
    public dbnsFlux
{

public:
    TypeName("rotatedFlux");
    
    // Constructors
    rotatedFlux() = delete;
    rotatedFlux(const rotatedFlux& ) = delete;
    
    //- Construct null
    rotatedFlux(const fvMesh&, const fluidThermo& thermo, const dictionary&);
    
    // Destructor - default
    virtual ~rotatedFlux() {}
    
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
        const scalar& RLeft,
        const scalar& RRight,
        const scalar& CvLeft,
        const scalar& CvRight,
        const vector& Sf,
        const scalar& magSf,
	const scalar& meshPhi
    ) const;

private:
    scalar epsilon_;
    autoPtr<dbnsFlux> diffusiveFlux_;
    autoPtr<dbnsFlux> lowDiffusionFlux_;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
