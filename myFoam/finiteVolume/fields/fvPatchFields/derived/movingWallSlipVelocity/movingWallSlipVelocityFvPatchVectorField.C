/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "movingWallSlipVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "symmTransformField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingWallSlipVelocityFvPatchVectorField::
movingWallSlipVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


Foam::movingWallSlipVelocityFvPatchVectorField::
movingWallSlipVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF)
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::movingWallSlipVelocityFvPatchVectorField::
movingWallSlipVelocityFvPatchVectorField
(
    const movingWallSlipVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::movingWallSlipVelocityFvPatchVectorField::
movingWallSlipVelocityFvPatchVectorField
(
    const movingWallSlipVelocityFvPatchVectorField& mwvpvf
)
:
    fixedValueFvPatchVectorField(mwvpvf)
{}


Foam::movingWallSlipVelocityFvPatchVectorField::
movingWallSlipVelocityFvPatchVectorField
(
    const movingWallSlipVelocityFvPatchVectorField& mwvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(mwvpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::movingWallSlipVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = internalField().mesh();

    if (mesh.moving())
    {
        const fvPatch& p = patch();

        const volVectorField& U =
            static_cast<const volVectorField&>(internalField());

        scalarField phip
        (
            p.patchField<surfaceScalarField, scalar>(fvc::meshPhi(U))
        );

        const vectorField n(p.nf());
        const scalarField& magSf = p.magSf();
        tmp<scalarField> Un = phip/(magSf + VSMALL);
	
	vectorField Urel(this->internalField() - Un()*n);

	vectorField::operator=
	  (
	   0.5*(Urel + transform(I - 2.0*sqr(n), Urel)) + Un()*n
 	   );
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::movingWallSlipVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
#if (OPENFOAM >= 1912)
    this->writeEntry("value", os);
#else
    writeEntry(os, "value", *this);
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        movingWallSlipVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
