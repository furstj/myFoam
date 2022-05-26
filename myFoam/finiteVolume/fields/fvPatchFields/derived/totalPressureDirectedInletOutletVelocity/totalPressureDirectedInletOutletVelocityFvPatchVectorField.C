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

#include "totalPressureDirectedInletOutletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::totalPressureDirectedInletOutletVelocityFvPatchVectorField::
totalPressureDirectedInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(p, iF),
    psiName_("none"),
    phiName_("phi"),
    gamma_(0.0),
    p0_(p.size(), 0.0),
    inletDir_(p.size())
{
    refValue() = *this;
    refGrad() = vector::zero;
    valueFraction() = 0.0;
}


Foam::totalPressureDirectedInletOutletVelocityFvPatchVectorField::
totalPressureDirectedInletOutletVelocityFvPatchVectorField
(
    const totalPressureDirectedInletOutletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchVectorField(ptf, p, iF, mapper),
    psiName_(ptf.psiName_),
    phiName_(ptf.phiName_),
    gamma_(ptf.gamma_),
    p0_(mapper(ptf.p0_)),
    inletDir_(mapper(ptf.inletDir_))
{}


Foam::totalPressureDirectedInletOutletVelocityFvPatchVectorField::
totalPressureDirectedInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchVectorField(p, iF),
    psiName_(dict.lookupOrDefault<word>("psi", "thermo:psi")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    gamma_(readScalar(dict.lookup("gamma"))),
    p0_("p0", dict, p.size()),
    inletDir_("inletDirection", dict, p.size())
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
    refValue() = *this;
    refGrad() = vector::zero;
    valueFraction() = 0.0;
}


Foam::totalPressureDirectedInletOutletVelocityFvPatchVectorField::
totalPressureDirectedInletOutletVelocityFvPatchVectorField
(
    const totalPressureDirectedInletOutletVelocityFvPatchVectorField& sfspvf,
    const DimensionedField<vector, volMesh>& iF
)
:
  mixedFvPatchVectorField(sfspvf, iF),
  psiName_(sfspvf.psiName_),
  phiName_(sfspvf.phiName_),
  gamma_(sfspvf.gamma_),
  p0_(sfspvf.p0_),
  inletDir_(sfspvf.inletDir_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::totalPressureDirectedInletOutletVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchVectorField::autoMap(m);
#if (OPENFOAM >= 1812)
    p0_.autoMap(m);
    inletDir_.autoMap(m);
#else
    m(p0_, p0_);
    m(inletDir_, inletDir_);
#endif
}


void Foam::totalPressureDirectedInletOutletVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    mixedFvPatchVectorField::rmap(ptf, addr);

    const totalPressureDirectedInletOutletVelocityFvPatchVectorField& tiptf =
        refCast<const totalPressureDirectedInletOutletVelocityFvPatchVectorField>(ptf);

    p0_.rmap(tiptf.p0_, addr);
    inletDir_.rmap(tiptf.inletDir_, addr);
}

void Foam::totalPressureDirectedInletOutletVelocityFvPatchVectorField::updateCoeffs()
{
  if (updated())
    {
      return;
    }

  const surfaceScalarField& phi =
    db().lookupObject<surfaceScalarField>(phiName_);
  
  const fvsPatchField<scalar>& phip =
    patch().patchField<surfaceScalarField, scalar>(phi);

  const fvPatchField<scalar>& psip =
    patch().lookupPatchField<volScalarField, scalar>(psiName_);

  const fvPatchField<scalar>& p =
    patch().lookupPatchField<volScalarField, scalar>("p");

  vectorField n = patch().nf();
  scalarField cosA = - n & inletDir_ / mag(inletDir_);
  
  scalar gM1ByG = (gamma_ - 1.0)/gamma_;

  scalarField sqrM = 2.0 / (gamma_ - 1) * ( pow(max(p0_/ p, 1.0), gM1ByG) - 1.0 );
 
  //refValue() = inletDir_ * sqrt( 2.0 * (pow(max(p0_ / p, 1.0), gM1ByG) - 1.0) / psip / gM1ByG);
  refValue() = inletDir_ * sqrt( min(sqrM, 1.0/sqr(cosA)) * gamma_ / psip );

  valueFraction() = 1.0 - pos(phip);

  mixedFvPatchVectorField::updateCoeffs();

}


void Foam::totalPressureDirectedInletOutletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    #if (OPENFOAM_PLUS>=1712 || OPENFOAM >= 1812)
    os.writeEntryIfDifferent<word>("psi", "thermo:psi", psiName_);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    #else
    writeEntryIfDifferent<word>(os, "psi", "thermo:psi", psiName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    #endif

    os.writeKeyword("gamma") << gamma_ << token::END_STATEMENT << nl;
#if (OPENFOAM >= 1812)
    p0_.writeEntry("p0", os);
    inletDir_.writeEntry("inletDirection", os);
    this->writeEntry("value", os);
#else
    writeEntry(os, "p0", p0_);
    writeEntry(os, "inletDirection", inletDir_);
    writeEntry(os, "value", *this);
#endif
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::totalPressureDirectedInletOutletVelocityFvPatchVectorField::operator=
(
 const fvPatchField<vector>& pvf
 )
{
  fvPatchField<vector>::operator=(inletDir_*(inletDir_ & pvf));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        totalPressureDirectedInletOutletVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
