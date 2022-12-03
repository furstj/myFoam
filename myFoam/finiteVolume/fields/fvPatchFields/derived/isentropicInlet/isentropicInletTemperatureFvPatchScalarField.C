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

#include "isentropicInletTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "isentropicInletPressureFvPatchScalarField.H"
#include "totalTemperatureFvPatchScalarField.H"          // <== TODO
#include "psiThermo.H"
#include "gasProperties.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isentropicInletTemperatureFvPatchScalarField::
isentropicInletTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    pName_("p"),
    T0_(p.size(), Zero)
{
}


Foam::isentropicInletTemperatureFvPatchScalarField::
isentropicInletTemperatureFvPatchScalarField
(
    const isentropicInletTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    pName_(ptf.pName_),
    T0_(ptf.T0_, mapper)
{}


Foam::isentropicInletTemperatureFvPatchScalarField::
isentropicInletTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    pName_(dict.getOrDefault<word>("p", "p")),
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


Foam::isentropicInletTemperatureFvPatchScalarField::
isentropicInletTemperatureFvPatchScalarField
(
    const isentropicInletTemperatureFvPatchScalarField& pivpvf
)
:
    fixedValueFvPatchScalarField(pivpvf),
    pName_(pivpvf.pName_),
    T0_(pivpvf.T0_)
{}


Foam::isentropicInletTemperatureFvPatchScalarField::
isentropicInletTemperatureFvPatchScalarField
(
    const isentropicInletTemperatureFvPatchScalarField& pivpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(pivpvf, iF),
    pName_(pivpvf.pName_),
    T0_(pivpvf.T0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::isentropicInletTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
#if (OPENFOAM >= 1812)
    T0_.autoMap(m);
#else
    m(T0_, T0_);
#endif    
}


void Foam::isentropicInletTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const isentropicInletTemperatureFvPatchScalarField& tiptf =
        refCast<const isentropicInletTemperatureFvPatchScalarField>(ptf);

    T0_.rmap(tiptf.T0_, addr);
}


void Foam::isentropicInletTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const auto& thermo =
        mesh.lookupObject<psiThermo>("thermophysicalProperties");

    autoPtr<gasProperties> gasProps(gasProperties::New(thermo));
    
    const fvPatchScalarField& pp = 
        patch().patchField<volScalarField, scalar>(
            db().lookupObjectRef<volScalarField>(pName_)
        );
    
    scalarField p0;
    if (pp.type() == "isentropicInletPressure")
    {
        p0 = refCast<const isentropicInletPressureFvPatchScalarField>(pp).p0();
    }
    else
    {
        FatalErrorIn("isentropicInletTemperatureFvPatchScalarField::updateCoeffs()") 
            << "the isentropicInletTemperature has to be combined with "
            << "isentropicInletPressure condition!"
            << abort(FatalError);
    }
    
    scalarField& pT = *this;

    forAll(pT, faceI) {
        scalar S = gasProps->S(p0[faceI], T0_[faceI]);
        scalar p = min(pp[faceI], p0[faceI]);
        scalar T = gasProps->TpS(p, S, pT[faceI]);
        pT[faceI] = T;
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::isentropicInletTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
#if (OPENFOAM_PLUS>=1712 || OPENFOAM >= 1812)
    os.writeEntryIfDifferent<word>("p", "p", pName_);
 #else
    writeEntryIfDifferent<word>(os, "p", "p", pName_);
 #endif
#if (OPENFOAM >= 1812)
    T0_.writeEntry("T0", os);
    this->writeEntry("value", os);
#else
    writeEntry(os, "T0", T0_);
    writeEntry(os, "value", *this);
#endif    

}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        isentropicInletTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
