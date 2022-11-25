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

#include "isentropicInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "isentropicInletPressureFvPatchScalarField.H"
#include "isentropicInletTemperatureFvPatchScalarField.H"
#include "psiThermo.H"
#include "gasProperties.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isentropicInletVelocityFvPatchVectorField::
isentropicInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(p, iF),
    pName_("p"),
    TName_("T"),
    inletDir_(p.size())
{
    refValue() = *this;
    refGrad() = Zero;
    valueFraction() = 0.0;
}


Foam::isentropicInletVelocityFvPatchVectorField::
isentropicInletVelocityFvPatchVectorField
(
    const isentropicInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchVectorField(ptf, p, iF, mapper),
    pName_(ptf.pName_),
    TName_(ptf.TName_),
    inletDir_(ptf.inletDir_, mapper)
{}


Foam::isentropicInletVelocityFvPatchVectorField::
isentropicInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchVectorField(p, iF),
    pName_(dict.getOrDefault<word>("p", "p")),
    TName_(dict.getOrDefault<word>("T", "T")),
    inletDir_("inletDirection", dict, p.size())
{
    patchType() = dict.getOrDefault<word>("patchType", word::null);
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
    refValue() = *this;
    refGrad() = Zero;
    valueFraction() = 0.0;
}


Foam::isentropicInletVelocityFvPatchVectorField::
isentropicInletVelocityFvPatchVectorField
(
    const isentropicInletVelocityFvPatchVectorField& pivpvf
)
:
    mixedFvPatchVectorField(pivpvf),
    pName_(pivpvf.pName_),
    TName_(pivpvf.TName_),
    inletDir_(pivpvf.inletDir_)
{}


Foam::isentropicInletVelocityFvPatchVectorField::
isentropicInletVelocityFvPatchVectorField
(
    const isentropicInletVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(pivpvf, iF),
    pName_(pivpvf.pName_),
    TName_(pivpvf.TName_),
    inletDir_(pivpvf.inletDir_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::isentropicInletVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchVectorField::autoMap(m);
#if (OPENFOAM >= 1812)
    inletDir_.autoMap(m);
#else
    m(inletDir_, inletDir_);
#endif    
}


void Foam::isentropicInletVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    mixedFvPatchVectorField::rmap(ptf, addr);

    const isentropicInletVelocityFvPatchVectorField& tiptf =
        refCast<const isentropicInletVelocityFvPatchVectorField>
        (ptf);

    inletDir_.rmap(tiptf.inletDir_, addr);
}


void Foam::isentropicInletVelocityFvPatchVectorField::updateCoeffs()
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
            db().lookupObject<volScalarField>(pName_)
        );
    const_cast<fvPatchScalarField&>(pp).updateCoeffs();
    
    const fvPatchScalarField& pT = 
        patch().patchField<volScalarField, scalar>(
            db().lookupObject<volScalarField>(TName_)
        );


    scalarField T0, p0;
    if (pT.type() == "isentropicInletTemperature" && pp.type() == "isentropicInletPressure")
    {
        T0 = refCast<const isentropicInletTemperatureFvPatchScalarField>(pT).T0();
        p0 = refCast<const isentropicInletPressureFvPatchScalarField>(pp).p0();
    }
    else
    {
        FatalErrorIn("isentropicInletVelocityFvPatchVectorField::updateCoeffs()") 
            << "the isentropicInletVelocity has to be combined with "
            << "isentropicInletTemperature and isentropicInletPressure conditions!"
            << abort(FatalError);
    }

    vectorField& refValue = this->refValue();
    scalarField& valFraction = this->valueFraction();

    forAll(refValue, faceI) {
        vector dir = inletDir_[faceI] / mag(inletDir_[faceI]);

        scalar S = gasProps->S(p0[faceI], T0[faceI]);
        scalar p = min(pp[faceI], p0[faceI]);

        scalar T = gasProps->TpS(p, S, pT[faceI]);

        scalar h  = gasProps->Hs(p, T);
        scalar H0 = gasProps->Hs(p0[faceI], T0[faceI]);

        scalar uMag = sqrt(2*max(H0 - h, 0.0));
        
        refValue[faceI] = dir*uMag;
        valFraction[faceI] = pos0(H0-h);
    }
    
    mixedFvPatchVectorField::updateCoeffs();
}


void Foam::isentropicInletVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
#if (OPENFOAM_PLUS>=1712 || OPENFOAM >= 1812)
    os.writeEntryIfDifferent<word>("p", "p", pName_);
    os.writeEntryIfDifferent<word>("T", "T", TName_);
 #else
    writeEntryIfDifferent<word>(os, "p", "p", pName_);
    writeEntryIfDifferent<word>(os, "T", "T", TName_);
 #endif
#if (OPENFOAM >= 1812)
    inletDir_.writeEntry("inletDirection", os);
    this->writeEntry("value", os);
#else
    writeEntry(os, "inletDirection", inletDir_);
    writeEntry(os, "value", *this);
#endif    

}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::isentropicInletVelocityFvPatchVectorField::operator=
(
    const fvPatchField<vector>& pvf
)
{
    fvPatchField<vector>::operator=
    (
        valueFraction()*(inletDir_*(inletDir_ & pvf))
      + (1 - valueFraction())*pvf
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        isentropicInletVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
