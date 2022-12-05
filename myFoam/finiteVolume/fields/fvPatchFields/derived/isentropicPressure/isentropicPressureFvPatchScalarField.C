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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isentropicPressureFvPatchScalarField::isentropicPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_("U"),
    phiName_("phi"),
    rhoName_("rho"),
    psiName_("none"),
    gamma_(0.0),
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
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    rhoName_(dict.getOrDefault<word>("rho", "rho")),
    psiName_(dict.getOrDefault<word>("psi", "none")),
    gamma_(psiName_ != "none" ? dict.get<scalar>("gamma") : 1),
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
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    psiName_(ptf.psiName_),
    gamma_(ptf.gamma_),
    p0_(ptf.p0_, mapper)
{}


Foam::isentropicPressureFvPatchScalarField::isentropicPressureFvPatchScalarField
(
    const isentropicPressureFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    psiName_(tppsf.psiName_),
    gamma_(tppsf.gamma_),
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
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    psiName_(tppsf.psiName_),
    gamma_(tppsf.gamma_),
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

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    if (internalField().dimensions() == dimPressure)
    {
        if (psiName_ == "none")
        {
            // Variable density and low-speed compressible flow

            const fvPatchField<scalar>& rho =
                patch().lookupPatchField<volScalarField, scalar>(rhoName_);

            operator==(p0p - 0.5*rho*(1.0 - pos0(phip))*magSqr(Up));
        }
        else
        {
            // High-speed compressible flow

            const fvPatchField<scalar>& psip =
                patch().lookupPatchField<volScalarField, scalar>(psiName_);

            if (gamma_ > 1)
            {
                scalar gM1ByG = (gamma_ - 1)/gamma_;

                operator==
                (
                    p0p
                   /pow
                    (
                        (1.0 + 0.5*psip*gM1ByG*(1.0 - pos0(phip))*magSqr(Up)),
                        1.0/gM1ByG
                    )
                );
            }
            else
            {
                operator==(p0p/(1.0 + 0.5*psip*(1.0 - pos0(phip))*magSqr(Up)));
            }
        }

    }
    else if (internalField().dimensions() == dimPressure/dimDensity)
    {
        // Incompressible flow
        operator==(p0p - 0.5*(1.0 - pos0(phip))*magSqr(Up));
    }
    else
    {
        FatalErrorInFunction
            << " Incorrect pressure dimensions " << internalField().dimensions()
            << nl
            << "    Should be " << dimPressure
            << " for compressible/variable density flow" << nl
            << "    or " << dimPressure/dimDensity
            << " for incompressible flow," << nl
            << "    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
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
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntry("rho", rhoName_);
    os.writeEntry("psi", psiName_);
    os.writeEntry("gamma", gamma_);
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
