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
    ausmPWPlusFlux

Description
    AUSMPW+ flux class, computes numerical fluxes based on AUSMPW+ scheme. Uses primitive 
    variables as input and gives back conservative numerical fluxes.

    KIM, Kyu Hong, KIM, Chongam and RHO, Oh Hyun. 
    Methods for the accurate computations of hypersonic flows. I. AUSMPW+ scheme. 
    Journal of Computational Physics. 2001. DOI 10.1006/jcph.2001.6873. 

Author
    Jiri Furstr  All rights reserved.


SourceFiles
    AUSMPWplusFlux.H

\*---------------------------------------------------------------------------*/

#ifndef AUSMPWplusFlux_H
#define AUSMPWplusFlux_H

#include "dbnsFlux.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class AUSMPWplusFlux Declaration
\*---------------------------------------------------------------------------*/

class AUSMPWplusFlux: public dbnsFlux
{
private:
    AUSMPWplusFlux();
    AUSMPWplusFlux(const AUSMPWplusFlux&);
  
public:
    TypeName("AUSMPWplusFlux");

    // Constructors

    //- Constructor 
    AUSMPWplusFlux(const fvMesh&, const fluidThermo& thermo, const dictionary&);

    // Destructor - default
    virtual ~AUSMPWplusFlux() {};

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
  
    scalar alpha_;
    scalar beta_;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
