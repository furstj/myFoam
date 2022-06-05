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

#include "totalPressureDirectedInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::totalPressureDirectedInletVelocityFvPatchVectorField::
totalPressureDirectedInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    psiName_("none"),
    gamma_(0.0),
    p0_(p.size(), 0.0),
    inletDir_(p.size())
{
}


Foam::totalPressureDirectedInletVelocityFvPatchVectorField::
totalPressureDirectedInletVelocityFvPatchVectorField
(
    const totalPressureDirectedInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    psiName_(ptf.psiName_),
    gamma_(ptf.gamma_),
    p0_(mapper(ptf.p0_)),
    inletDir_(mapper(ptf.inletDir_))
{}


Foam::totalPressureDirectedInletVelocityFvPatchVectorField::
totalPressureDirectedInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    psiName_(dict.lookupOrDefault<word>("psi", "thermo:psi")),
    gamma_(readScalar(dict.lookup("gamma"))),
    p0_("p0", dict, p.size()),
    inletDir_("inletDirection", dict, p.size())
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}

Foam::totalPressureDirectedInletVelocityFvPatchVectorField::
totalPressureDirectedInletVelocityFvPatchVectorField
(
    const totalPressureDirectedInletVelocityFvPatchVectorField& sfspvf
)
  :
  fixedValueFvPatchVectorField(sfspvf),
  psiName_(sfspvf.psiName_),
  gamma_(sfspvf.gamma_),
  p0_(sfspvf.p0_),
  inletDir_(sfspvf.inletDir_)
{}


Foam::totalPressureDirectedInletVelocityFvPatchVectorField::
totalPressureDirectedInletVelocityFvPatchVectorField
(
    const totalPressureDirectedInletVelocityFvPatchVectorField& sfspvf,
    const DimensionedField<vector, volMesh>& iF
)
:
  fixedValueFvPatchVectorField(sfspvf, iF),
  psiName_(sfspvf.psiName_),
  gamma_(sfspvf.gamma_),
  p0_(sfspvf.p0_),
  inletDir_(sfspvf.inletDir_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::totalPressureDirectedInletVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
#if (OPENFOAM >= 1812)
    p0_.autoMap(m);
    inletDir_.autoMap(m);
#else
    m(p0_, p0_);
    m(inletDir_, inletDir_);
#endif
}


void Foam::totalPressureDirectedInletVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const totalPressureDirectedInletVelocityFvPatchVectorField& tiptf =
        refCast<const totalPressureDirectedInletVelocityFvPatchVectorField>(ptf);

    p0_.rmap(tiptf.p0_, addr);
    inletDir_.rmap(tiptf.inletDir_, addr);
}

void Foam::totalPressureDirectedInletVelocityFvPatchVectorField::updateCoeffs()
{
  if (updated())
    {
      return;
    }

  const fvPatchField<scalar>& psip =
    patch().lookupPatchField<volScalarField, scalar>(psiName_);

  const fvPatchField<scalar>& p =
    patch().lookupPatchField<volScalarField, scalar>("p");

  vectorField n = patch().nf();
  scalarField cosA = - n & inletDir_ / mag(inletDir_);
  
  scalar gM1ByG = (gamma_ - 1.0)/gamma_;

  scalarField sqrM = 2.0 / (gamma_ - 1) * ( pow(max(p0_/ p, 1.0), gM1ByG) - 1.0 );
 
  operator==(inletDir_ * sqrt( min(sqrM, 1.0/sqr(cosA)) * gamma_ / psip ));

  fixedValueFvPatchVectorField::updateCoeffs();

}


void Foam::totalPressureDirectedInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    #if (OPENFOAM_PLUS>=1712 || OPENFOAM >= 1812)
    os.writeEntryIfDifferent<word>("psi", "thermo:psi", psiName_);
    #else
    writeEntryIfDifferent<word>(os, "psi", "thermo:psi", psiName_);
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

void Foam::totalPressureDirectedInletVelocityFvPatchVectorField::operator=
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
        totalPressureDirectedInletVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
