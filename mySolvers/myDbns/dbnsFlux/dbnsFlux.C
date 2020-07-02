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
    dbnsFlux.C

\*---------------------------------------------------------------------------*/

#include "dbnsFlux.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  defineTypeNameAndDebug(dbnsFlux, 0);
  defineRunTimeSelectionTable(dbnsFlux,dictionary);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::dbnsFlux> Foam::dbnsFlux::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word key
)
{
    const word fluxName( dict.lookup(key) );

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(fluxName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "dbnsFlux::New(mesh, dict, key)"
        )   << "Unknown dbnsFlux type "
            << fluxName << nl << nl
            << "Valid dbnsFlux types are :" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<dbnsFlux>(cstrIter()(mesh, dict));
}


// ************************************************************************* //
