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
#include "patchTransformedInterpolation.H"
#include "pointFields.H"
#include "symmTensor2D.H"
#include "tensor2D.H"
#include "syncTools.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchTransformedInterpolation, 0);

    addToRunTimeSelectionTable
    (
        motionInterpolation,
        patchTransformedInterpolation,
        Istream
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::patchTransformedInterpolation::getPatches
(
    Istream& entry
) const
{
    wordList patchNames(entry);

    labelList patches(patchNames.size(), -1);

    forAll(patchNames, patchI)
    {
        patches[patchI] =
            mesh().boundaryMesh().findPatchID
            (
                patchNames[patchI]
            );

        if (patches[patchI] == -1)
        {
            FatalErrorInFunction
                << "patch \"" << patchNames[patchI]
                << "\" not found" << exit(FatalError);
        }
    }

    return patches;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchTransformedInterpolation::patchTransformedInterpolation
(
    const fvMesh& mesh,
    Istream& entry
)
:
    motionInterpolation(mesh, entry),
    patches_(getPatches(entry))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchTransformedInterpolation::~patchTransformedInterpolation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::patchTransformedInterpolation::interpolate
(
    const volScalarField&,
    pointScalarField&
) const
{
    NotImplemented;
}


void Foam::patchTransformedInterpolation::interpolate
(
    const volVectorField& cellDisplacement,
    pointVectorField& pointDisplacement
) const
{
    const pointField& points(mesh().points());
    const label nPoints(points.size());

    volPointInterpolation::New(mesh()).interpolate
    (
        cellDisplacement,
        pointDisplacement
    );

    pointDisplacement.correctBoundaryConditions();

    vectorField pointRotation(nPoints, Zero);
    scalarField pointExpansion(nPoints, scalar(0));

    labelList pointDisplacementNSum(nPoints, 0);
    vectorField pointDisplacementSum(nPoints, Zero);

    forAll(patches_, patchI)
    {
        const polyPatch& patch(mesh().boundaryMesh()[patches_[patchI]]);

        forAll(patch, pFaceI)
        {
            const face& f(patch[pFaceI]);

            const label cellI(patch.faceCells()[pFaceI]);
            const cell& c(mesh().cells()[cellI]);
            const labelList cPoints(c.labels(mesh().faces()));

            // Consider movement around the face centre
            const point& xOrigin(patch.faceCentres()[pFaceI]);

            // Mean translation
            const vector uMean(f.average(points, pointDisplacement));

            // Calculate rotation and expansion for each point
            forAll(f, fPointI)
            {
                const label pointI(f[fPointI]);
                const vector& x(points[pointI]);
                const vector r(x - xOrigin);
                const vector u(pointDisplacement[pointI] - uMean);

                pointRotation[pointI] = 2*(r ^ u)/magSqr(r);
                pointExpansion[pointI] = (r & u)/magSqr(r);
            }

            // Mean rotation and expansion
            const vector omegaMean(f.average(points, pointRotation));
            const scalar divMean(f.average(points, pointExpansion));

            // Apply mean solid body motion to all cell points
            forAll(cPoints, cPointI)
            {
                const label pointI(cPoints[cPointI]);
                const vector& x(points[pointI]);
                const vector r(x - xOrigin);

                pointDisplacementNSum[pointI] += 1;
                pointDisplacementSum[pointI] +=
                    uMean + (omegaMean ^ r) + (divMean*r);
            }
        }
    }

    syncTools::syncPointList
    (
        mesh(),
        pointDisplacementNSum,
        plusEqOp<label>(),
        label(0)
    );

    syncTools::syncPointList
    (
        mesh(),
        pointDisplacementSum,
        plusEqOp<vector>(),
        vector::zero
    );

    forAll(points, pointI)
    {
        if (pointDisplacementNSum[pointI])
        {
            pointDisplacement[pointI] =
                pointDisplacementSum[pointI]/pointDisplacementNSum[pointI];
        }
    }

    // Correct the faces
    pointDisplacement.correctBoundaryConditions();
}


// ************************************************************************* //
