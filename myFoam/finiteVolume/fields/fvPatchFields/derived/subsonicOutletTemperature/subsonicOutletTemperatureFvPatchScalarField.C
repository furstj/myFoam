/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
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

#include "subsonicOutletTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "psiThermo.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::subsonicOutletTemperatureFvPatchScalarField::
subsonicOutletTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    pName_("p"),
    rhoName_("rho"),
    psiName_("thermo:psi")
{
    refValue() = *this;
    refGrad() = Zero;
    valueFraction() = 0.0;
}


Foam::subsonicOutletTemperatureFvPatchScalarField::
subsonicOutletTemperatureFvPatchScalarField
(
    const subsonicOutletTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    pName_(ptf.pName_),
    rhoName_(ptf.rhoName_),
    psiName_(ptf.psiName_)
{}


Foam::subsonicOutletTemperatureFvPatchScalarField::
subsonicOutletTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    pName_(dict.lookupOrDefault<word>("p", "p")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    psiName_(dict.lookupOrDefault<word>("psi", "thermo:psi"))
{
    refValue() = *this;
    refGrad() = Zero;
    valueFraction() = 0.0;

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
            (
                scalarField("value", dict, p.size())
            );
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
    }
}


Foam::subsonicOutletTemperatureFvPatchScalarField::
subsonicOutletTemperatureFvPatchScalarField
(
    const subsonicOutletTemperatureFvPatchScalarField& sfspvf
)
:
    mixedFvPatchScalarField(sfspvf),
    pName_(sfspvf.pName_),
    rhoName_(sfspvf.rhoName_),
    psiName_(sfspvf.psiName_)
{}


Foam::subsonicOutletTemperatureFvPatchScalarField::
subsonicOutletTemperatureFvPatchScalarField
(
    const subsonicOutletTemperatureFvPatchScalarField& sfspvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(sfspvf, iF),
    pName_(sfspvf.pName_),
    rhoName_(sfspvf.rhoName_),
    psiName_(sfspvf.psiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::subsonicOutletTemperatureFvPatchScalarField::updateCoeffs()
{
    if (!size() || updated())
    {
        return;
    }
    
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const fluidThermo& thermo =
        mesh.lookupObject<fluidThermo>("thermophysicalProperties");

    tmp< volScalarField > gamma = thermo.gamma();
    const fvPatchField<scalar>&  pgamma =
        gamma->boundaryField()[patch().index()];

    tmp< volScalarField > Cp = thermo.Cp();
    const fvPatchField<scalar>&  pCp =
        Cp->boundaryField()[patch().index()];

    tmp< volScalarField > Cv = thermo.Cv();
    const fvPatchField<scalar>&  pCv =
        Cv->boundaryField()[patch().index()];

    const fvPatchField<scalar>& ppsi =
        thermo.psi().boundaryField()[patch().index()];

    const fvPatchField<scalar>& pp =
        patch().lookupPatchField<volScalarField, scalar>(pName_);

    const fvPatchField<scalar>& prho =
        patch().lookupPatchField<volScalarField, scalar>(rhoName_);

    valueFraction() = 0;
    refGrad() = Zero;

    const scalarField pi(pp.patchInternalField());
    const scalarField rhoi(prho.patchInternalField());
    const scalarField ciSqr(pgamma.patchInternalField()/ppsi.patchInternalField());

    scalarField rhop(rhoi + (pp - pi)/ciSqr);

    this->operator==(pp/rhop/(pCp - pCv));
}


void Foam::subsonicOutletTemperatureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    #if (OPENFOAM_PLUS>=1712 || OPENFOAM >=1812)
    os.writeEntryIfDifferent<word>("p", "p", pName_);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
    os.writeEntryIfDifferent<word>("psi", "thermo:psi", psiName_);
    #else
    writeEntryIfDifferent<word>(os, "p", "p", pName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntryIfDifferent<word>(os, "psi", "thermo:psi", psiName_);
    #endif
    #if (OPENFOAM >= 1812)
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
        fvPatchScalarField,
        subsonicOutletTemperatureFvPatchScalarField
    );
}



// ************************************************************************* //
