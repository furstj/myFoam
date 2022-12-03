/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "isentropicInletPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isentropicInletPressureFvPatchScalarField::
isentropicInletPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    p0_(p.size(), Zero)
{
    refValue() = *this;
    refGrad() = Zero;
    valueFraction() = 0.0;
}


Foam::isentropicInletPressureFvPatchScalarField::
isentropicInletPressureFvPatchScalarField
(
    const isentropicInletPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    p0_(ptf.p0_, mapper)
{}


Foam::isentropicInletPressureFvPatchScalarField::
isentropicInletPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    p0_("p0", dict, p.size())
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(p0_);
    }

    refGrad() = Zero;
    valueFraction() = 0.0;
}


Foam::isentropicInletPressureFvPatchScalarField::
isentropicInletPressureFvPatchScalarField
(
    const isentropicInletPressureFvPatchScalarField& pivpvf
)
:
    mixedFvPatchScalarField(pivpvf),
    p0_(pivpvf.p0_)
{}


Foam::isentropicInletPressureFvPatchScalarField::
isentropicInletPressureFvPatchScalarField
(
    const isentropicInletPressureFvPatchScalarField& pivpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(pivpvf, iF),
    p0_(pivpvf.p0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::isentropicInletPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
#if (OPENFOAM >= 1812)
    p0_.autoMap(m);
#else
    m(p0_, p0_);
#endif    
}


void Foam::isentropicInletPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const isentropicInletPressureFvPatchScalarField& tiptf =
        refCast<const isentropicInletPressureFvPatchScalarField>
        (ptf);

    p0_.rmap(tiptf.p0_, addr);
}


void Foam::isentropicInletPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    mixedFvPatchScalarField::updateCoeffs();    
}


void Foam::isentropicInletPressureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
#if (OPENFOAM >= 1812)
    p0_.writeEntry("p0", os);
    this->writeEntry("value", os);
#else
    writeEntry(os, "p0", p0_);
    writeEntry(os, "value", *this);
#endif    

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        isentropicInletPressureFvPatchScalarField
    );
}

// ************************************************************************* //
