/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: Open Source CFD
   \\    /   O peration     | 
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  | 
-------------------------------------------------------------------------------
License
    This file isn't part of foam-extend nor OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    dbnsFlux

Description
    Basic class for of inviscid numerical fluxes.

Author
    Jiri Furst

SourceFiles
    dbnsFlux.H

\*---------------------------------------------------------------------------*/

#ifndef dbnsFlux_H
#define dbnsFlux_H

#include "fvCFD.H"
#include "gasProperties.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class roeFlux Declaration
\*---------------------------------------------------------------------------*/

class dbnsFlux
{
private:
    dbnsFlux(const dbnsFlux&);
    void operator=(const dbnsFlux&);

protected:
    autoPtr<gasProperties> gasProps_;
    
public:
    TypeName("dbnsFlux");

    declareRunTimeSelectionTable
    (
        autoPtr,
        dbnsFlux,
        dictionary,
        (const fvMesh& mesh, const fluidThermo& thermo, const dictionary& dict),
        (mesh, thermo, dict)
    );

    // Constructors
    dbnsFlux() = delete;

    dbnsFlux(const fluidThermo& thermo):
        gasProps_(gasProperties::New(thermo))
    {};
    
    //- Return a reference to the selected dbnsFlux
    static autoPtr<dbnsFlux> New(
        const fvMesh& mesh,
        const fluidThermo& thermo,
        const dictionary& dict,
	const word key="dbnsFlux"
    );

    // Destructor - default
    virtual ~dbnsFlux() {};

    // Member Functions

        // Update fluxes based on current left and right states
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
        ) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
