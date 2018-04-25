/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"
#include "motionInterpolation.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(motionInterpolation, 0);

    defineRunTimeSelectionTable(motionInterpolation, Istream);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::motionInterpolation::motionInterpolation
(
    const fvMesh& mesh
)
:
    mesh_(mesh)
{}


Foam::motionInterpolation::motionInterpolation
(
    const fvMesh& mesh,
    Istream& entry
)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::motionInterpolation>
Foam::motionInterpolation::New(const fvMesh& mesh)
{
    return autoPtr<motionInterpolation>(new motionInterpolation(mesh));
}


Foam::autoPtr<Foam::motionInterpolation>
Foam::motionInterpolation::New(const fvMesh& mesh, Istream& entry)
{
    const word type(entry);

    Info<< "Selecting motion interpolation: " << type << endl;

    IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(type);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown interpolation type "
            << type << nl << nl
            << "Valid interpolation types are :" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<motionInterpolation>(cstrIter()(mesh, entry));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::motionInterpolation::~motionInterpolation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::motionInterpolation::interpolate
(
    const volScalarField& cellDisplacement,
    pointScalarField& pointDisplacement
) const
{
    volPointInterpolation::New(mesh()).interpolate
    (
        cellDisplacement,
        pointDisplacement
    );
}


void Foam::motionInterpolation::interpolate
(
    const volVectorField& cellDisplacement,
    pointVectorField& pointDisplacement
) const
{
    volPointInterpolation::New(mesh()).interpolate
    (
        cellDisplacement,
        pointDisplacement
    );
}


// ************************************************************************* //
