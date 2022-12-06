/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "isentropicPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fluidThermo.H"
#include "gasProperties.H"
#include "isentropicTemperatureFvPatchScalarField.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isentropicPressureFvPatchScalarField::isentropicPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),  
    UName_("U"),
    TName_("T"),
    phiName_("phi"),
    p0_(p.size(), Zero)
{}


Foam::isentropicPressureFvPatchScalarField::isentropicPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    UName_(dict.getOrDefault<word>("U", "U")),
    TName_(dict.getOrDefault<word>("T", "T")),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
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
}


Foam::isentropicPressureFvPatchScalarField::isentropicPressureFvPatchScalarField
(
    const isentropicPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    TName_(ptf.TName_),
    phiName_(ptf.phiName_),
    p0_(ptf.p0_, mapper)
{}


Foam::isentropicPressureFvPatchScalarField::isentropicPressureFvPatchScalarField
(
    const isentropicPressureFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    UName_(tppsf.UName_),
    TName_(tppsf.TName_),
    phiName_(tppsf.phiName_),
    p0_(tppsf.p0_)
{}


Foam::isentropicPressureFvPatchScalarField::isentropicPressureFvPatchScalarField
(
    const isentropicPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    UName_(tppsf.UName_),
    TName_(tppsf.TName_),
    phiName_(tppsf.phiName_),
    p0_(tppsf.p0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::isentropicPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    p0_.autoMap(m);
}


void Foam::isentropicPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const isentropicPressureFvPatchScalarField& tiptf =
        refCast<const isentropicPressureFvPatchScalarField>(ptf);

    p0_.rmap(tiptf.p0_, addr);
}


void Foam::isentropicPressureFvPatchScalarField::updateCoeffs
(
    const scalarField& p0p,
    const vectorField& Up
)
{
    if (updated())
    {
        return;
    }

    const fvPatchScalarField& Tp =
        patch().lookupPatchField<volScalarField, scalar>(TName_);

    const fvsPatchScalarField& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const fluidThermo& thermo =
        db().lookupObject<fluidThermo>(fluidThermo::dictName);

    autoPtr<gasProperties> gasProps(gasProperties::New(thermo));
    
    scalarField& pp = *this;

    const scalarField& T0 = refCast<const isentropicTemperatureFvPatchScalarField>(Tp).T0();
    
    forAll(Tp, faceI)
    {
        if (phip[faceI] < 0)
        {
            scalar S  = gasProps->S(p0p[faceI], T0[faceI]);
            scalar Hs = gasProps->Hs(p0p[faceI], T0[faceI]) - 0.5*magSqr(Up[faceI]);
            pp[faceI] = gasProps->pHS(Hs, S, pp[faceI]);
        }
        else
        {
            pp[faceI] = p0_[faceI];
        };
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::isentropicPressureFvPatchScalarField::updateCoeffs()
{
    updateCoeffs
    (
        p0(),
        patch().lookupPatchField<volVectorField, vector>(UName())
    );
}


void Foam::isentropicPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntryIfDifferent<word>("U", "U", UName_);
    os.writeEntryIfDifferent<word>("T", "T", TName_);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    p0_.writeEntry("p0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        isentropicPressureFvPatchScalarField
    );
}

// ************************************************************************* //
