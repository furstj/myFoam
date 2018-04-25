/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

#include "fixedValuePointPatchField.H"
#include "PointData.H"
#include "PointEdgeWave.H"
#include "volPointInterpolation.H"
#include "zeroGradientPointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template <class Type>
void Foam::patchCorrectedInterpolation::interpolateType
(
    const GeometricField<Type, fvPatchField, volMesh>& cellDisplacement,
    GeometricField<Type, pointPatchField, pointMesh>& pointDisplacement
) const
{
    // Create an uncorrected field
    GeometricField<Type, pointPatchField, pointMesh>
        pointUncorrectedDisplacement
        (
            IOobject
            (
                "pointUncorrectedDisplacement",
                mesh().time().timeName(),
                mesh()
            ),
            pointDisplacement.mesh(),
            pointDisplacement.dimensions(),
            fixedValuePointPatchField<Type>::typeName
        );

    // Interpolate to the uncorrected field, overwriting the fixed value
    // boundary conditions

    /*********************************************************************** 
        DISABLED due to bug in OF41
     *********************************************************************** 
    pointUncorrectedDisplacement ==
        volPointInterpolation::New(mesh()).interpolate
        (
            cellDisplacement,
            wordList
            (
                pointUncorrectedDisplacement.boundaryField().size(),
                zeroGradientPointPatchField<Type>::typeName
            )
        );

     *********************************************************************** 
        MANUAL FIX TO PREVIOUS BUG:
    */
    
    const pointMesh& pm = pointMesh::New(cellDisplacement.mesh());

    tmp<GeometricField<Type, pointPatchField, pointMesh>> tpf
        (
            new GeometricField<Type, pointPatchField, pointMesh>
            (
                IOobject
                (
                    "volPointInterpolate(" + cellDisplacement.name() + ')',
                    cellDisplacement.instance(),
                    pm.thisDb()
                ),
                pm,
                cellDisplacement.dimensions(),
                wordList
                (
                    pointUncorrectedDisplacement.boundaryField().size(),
                    zeroGradientPointPatchField<Type>::typeName
                )
            )
        );
    
    volPointInterpolation::New(mesh()).interpolateInternalField(cellDisplacement, tpf.ref());
    volPointInterpolation::New(mesh()).interpolateBoundaryField(cellDisplacement, tpf.ref(), true);

    pointUncorrectedDisplacement == tpf;
    /* END OF FIX */

    // Set the point displacement to the uncorrected result everywhere except
    // for on the boundary
    pointDisplacement.primitiveFieldRef() =
        pointUncorrectedDisplacement.primitiveField();
    pointDisplacement.correctBoundaryConditions();

    // Set the residual displacement as the difference between the boundary
    // specification and the uncorrected solution
    // (this reuses the uncorrected displacement field as the residual)
    pointUncorrectedDisplacement ==
        pointDisplacement - pointUncorrectedDisplacement;

    // Interpolate the residual from the boundary into the field
    interpolateDataFromPatchGroups(pointUncorrectedDisplacement);

    // Add the residual to the point displacement and correct the boundary
    pointDisplacement += pointUncorrectedDisplacement;
    pointDisplacement.correctBoundaryConditions();
}


template <class Type>
void Foam::patchCorrectedInterpolation::interpolateDataFromPatchGroups
(
    GeometricField<Type, pointPatchField, pointMesh>& data
) const
{
    // Initialise
    pointScalarField weight
    (
        IOobject
        (
            "weight",
            mesh().time().timeName(),
            mesh()
        ),
        data.mesh(),
        dimensionedScalar("zero", dimless, 0),
        zeroGradientPointPatchField<scalar>::typeName
    );
    data = dimensioned<Type>("zero", data.dimensions(), Type(Zero));

    forAll(patchGroups_, patchGroupI)
    {
        // Distance and data for the current group
        pointScalarField patchDistance
        (
            IOobject
            (
                "patchDistance",
                mesh().time().timeName(),
                mesh()
            ),
            data.mesh(),
            dimensionedScalar("zero", data.dimensions(), 0),
            zeroGradientPointPatchField<scalar>::typeName
        );
        GeometricField<Type, pointPatchField, pointMesh> patchData(data);

        // Wave the data through the mesh
        propagateDataFromPatchGroup
        (
            patchGroupI,
            patchDistance,
            patchData
        );

        // Calculate the weight and add to weighted sum
        const scalarField patchWeight
        (
            1/max(sqr(patchDistance.primitiveField()), SMALL)
        );
        data.primitiveFieldRef() += patchWeight*patchData.primitiveField();
        weight.primitiveFieldRef() += patchWeight;
    }

    // Complete the average
    data /= weight;
}


template <class Type>
void Foam::patchCorrectedInterpolation::propagateDataFromPatchGroup
(
    const label patchGroupI,
    pointScalarField& distance,
    GeometricField<Type, pointPatchField, pointMesh>& data
) const
{
    const labelList& patchGroup(patchGroups_[patchGroupI]);

    // Get the size of the seed info
    label nSeedInfo(0);
    forAll(patchGroup, patchGroupI)
    {
        const label patchI(patchGroup[patchGroupI]);

        nSeedInfo += data.boundaryField()[patchI].size();
    }

    // Generate the seed labels and info
    labelList seedLabels(nSeedInfo);
    List<PointData<Type>> seedInfo(nSeedInfo);
    nSeedInfo = 0;
    forAll(patchGroup, patchGroupI)
    {
        const label patchI(patchGroup[patchGroupI]);

        pointPatchField<Type>& patchDataField(data.boundaryFieldRef()[patchI]);

        patchDataField.updateCoeffs();

        const pointPatch& patch(patchDataField.patch());
        const Field<Type> patchData(patchDataField.patchInternalField());

        forAll(patch.meshPoints(), patchPointI)
        {
            const label pointI(patch.meshPoints()[patchPointI]);

            seedLabels[nSeedInfo] = pointI;

            seedInfo[nSeedInfo] =
                PointData<Type>
                (
                    mesh().points()[pointI],
                    0,
                    patchData[patchPointI]
                );

            nSeedInfo ++;
        }
    }

    // Wave the data through the mesh
    List<PointData<Type>> allPointInfo(mesh().nPoints());
    List<PointData<Type>> allEdgeInfo(mesh().nEdges());
    PointEdgeWave<PointData<Type>>
    (
        mesh(),
        seedLabels,
        seedInfo,
        allPointInfo,
        allEdgeInfo,
        mesh().globalData().nTotalPoints()
    );

    // Copy result into the fields
    forAll(allPointInfo, pointI)
    {
        distance[pointI] = sqrt(allPointInfo[pointI].distSqr());
        data[pointI] = allPointInfo[pointI].data();
    }
}

// ************************************************************************* //
