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

#include "uniformIsentropicPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fluidThermo.H"
#include "gasProperties.H"
#include "isentropicTemperatureFvPatchScalarField.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniformIsentropicPressureFvPatchScalarField::uniformIsentropicPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),  
    UName_("U"),
    TName_("T"),
    phiName_("phi"),
    p0_()
{}


Foam::uniformIsentropicPressureFvPatchScalarField::uniformIsentropicPressureFvPatchScalarField
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
    p0_(Function1<scalar>::New("p0", dict, &db()))
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
        const scalar t = this->db().time().timeOutputValue();
        fvPatchField<scalar>::operator=(p0_->value(t));
    }
}


Foam::uniformIsentropicPressureFvPatchScalarField::uniformIsentropicPressureFvPatchScalarField
(
    const uniformIsentropicPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    TName_(ptf.TName_),
    phiName_(ptf.phiName_),
    p0_(ptf.p0_.clone())
{
    patchType() = ptf.patchType();
    
    // This is not ideal but avoids problems with the creation of patch faces
    const scalar t = this->db().time().timeOutputValue();
    fvPatchScalarField::operator==(p0_->value(t));
}   // Set the patch pressure to the current total pressure



Foam::uniformIsentropicPressureFvPatchScalarField::uniformIsentropicPressureFvPatchScalarField
(
    const uniformIsentropicPressureFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    UName_(tppsf.UName_),
    TName_(tppsf.TName_),
    phiName_(tppsf.phiName_),
    p0_(tppsf.p0_.clone())
{}


Foam::uniformIsentropicPressureFvPatchScalarField::uniformIsentropicPressureFvPatchScalarField
(
    const uniformIsentropicPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    UName_(tppsf.UName_),
    TName_(tppsf.TName_),
    phiName_(tppsf.phiName_),
    p0_(tppsf.p0_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::uniformIsentropicPressureFvPatchScalarField::updateCoeffs
(
    const vectorField& Up
)
{
    if (updated())
    {
        return;
    }

    scalar p0 = p0_->value(this->db().time().timeOutputValue());
    
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
            scalar S  = gasProps->S(p0, T0[faceI]);
            scalar Hs = gasProps->Hs(p0, T0[faceI]) - 0.5*magSqr(Up[faceI]);
            pp[faceI] = gasProps->pHS(Hs, S, pp[faceI]);
        }
        else
        {
            pp[faceI] = p0;
        };
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::uniformIsentropicPressureFvPatchScalarField::updateCoeffs()
{
    updateCoeffs
    (
        patch().lookupPatchField<volVectorField, vector>(UName())
    );
}


void Foam::uniformIsentropicPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntryIfDifferent<word>("U", "U", UName_);
    os.writeEntryIfDifferent<word>("T", "T", TName_);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    p0_->writeData(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        uniformIsentropicPressureFvPatchScalarField
    );
}

// ************************************************************************* //
