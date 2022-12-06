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

#include "isentropicTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fluidThermo.H"
#include "gasProperties.H"
#include "isentropicPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isentropicTemperatureFvPatchScalarField::isentropicTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    pName_("p"),
    UName_("U"),
    phiName_("phi"),
    T0_(p.size(), Zero)
{}


Foam::isentropicTemperatureFvPatchScalarField::isentropicTemperatureFvPatchScalarField
(
    const isentropicTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    pName_(ptf.pName_),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    T0_(ptf.T0_, mapper)
{}


Foam::isentropicTemperatureFvPatchScalarField::isentropicTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    pName_(dict.getOrDefault<word>("p", "p")),
    UName_(dict.getOrDefault<word>("U", "U")),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    T0_("T0", dict, p.size())
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
        fvPatchField<scalar>::operator=(T0_);
    }
}


Foam::isentropicTemperatureFvPatchScalarField::isentropicTemperatureFvPatchScalarField
(
    const isentropicTemperatureFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    pName_(tppsf.pName_),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    T0_(tppsf.T0_)
{}


Foam::isentropicTemperatureFvPatchScalarField::isentropicTemperatureFvPatchScalarField
(
    const isentropicTemperatureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    pName_(tppsf.pName_),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    T0_(tppsf.T0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::isentropicTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    T0_.autoMap(m);
}


void Foam::isentropicTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const isentropicTemperatureFvPatchScalarField& tiptf =
        refCast<const isentropicTemperatureFvPatchScalarField>(ptf);

    T0_.rmap(tiptf.T0_, addr);
}


void Foam::isentropicTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchScalarField& pp =
        patch().lookupPatchField<volScalarField, scalar>(pName_);

    const fvPatchVectorField& Up =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    const fvsPatchScalarField& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const fluidThermo& thermo =
        db().lookupObject<fluidThermo>(fluidThermo::dictName);

    autoPtr<gasProperties> gasProps(gasProperties::New(thermo));
    
    scalarField& Tp = *this;

    const scalarField& p0 = refCast<const isentropicPressureFvPatchScalarField>(pp).p0();
    
    forAll(Tp, faceI)
    {
        if (phip[faceI] < 0)
        {
            scalar S  = gasProps->S(p0[faceI], T0_[faceI]);
            scalar Hs = gasProps->Hs(p0[faceI], T0_[faceI]) - 0.5*magSqr(Up[faceI]);
            scalar p = gasProps->pHS(Hs, S, pp[faceI]);
            Tp[faceI] = gasProps->TpS(p, S, Tp[faceI]);
        }
        else
        {
            Tp[faceI] = T0_[faceI];
        };
    }
    
    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::isentropicTemperatureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntryIfDifferent<word>("p", "p", pName_);
    os.writeEntryIfDifferent<word>("U", "U", UName_);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    T0_.writeEntry("T0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        isentropicTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
