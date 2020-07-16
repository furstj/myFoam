/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Author
    Sebastian Saegeler  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "mixedSubsonicSupersonicOutletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

#include "basicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mixedSubsonicSupersonicOutletFvPatchScalarField::
mixedSubsonicSupersonicOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    TName_("T"),
    fixedValue_(p.size(), 0.0),
    UName_("U")
{
    refValue() = *this;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


mixedSubsonicSupersonicOutletFvPatchScalarField::
mixedSubsonicSupersonicOutletFvPatchScalarField
(
    const mixedSubsonicSupersonicOutletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    fixedValue_(mapper(ptf.fixedValue_)),
    UName_(ptf.UName_)
{}


mixedSubsonicSupersonicOutletFvPatchScalarField::
mixedSubsonicSupersonicOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    TName_("T"),
    fixedValue_("fixedValue", dict, p.size()),
    UName_("U")
{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    refValue() = *this;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


mixedSubsonicSupersonicOutletFvPatchScalarField::
mixedSubsonicSupersonicOutletFvPatchScalarField
(
    const mixedSubsonicSupersonicOutletFvPatchScalarField& pivpvf
)
:
    mixedFvPatchScalarField(pivpvf),
    TName_(pivpvf.TName_),
    fixedValue_(pivpvf.fixedValue_),
    UName_(pivpvf.UName_)
{}


mixedSubsonicSupersonicOutletFvPatchScalarField::
mixedSubsonicSupersonicOutletFvPatchScalarField
(
    const mixedSubsonicSupersonicOutletFvPatchScalarField& pivpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(pivpvf, iF),
    TName_(pivpvf.TName_),
    fixedValue_(pivpvf.fixedValue_),
    UName_(pivpvf.UName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void mixedSubsonicSupersonicOutletFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);

#if (OPENFOAM >= 1912)
    fixedValue_.autoMap(m);
#else
    m(fixedValue_, fixedValue_);
#endif
}


void mixedSubsonicSupersonicOutletFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const mixedSubsonicSupersonicOutletFvPatchScalarField& tiptf =
        refCast<const mixedSubsonicSupersonicOutletFvPatchScalarField>(ptf);

    fixedValue_.rmap(tiptf.fixedValue_, addr);
}


void mixedSubsonicSupersonicOutletFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
        const fvPatchField<scalar>& Tp =
            patch().lookupPatchField<volScalarField, scalar>(TName_);
	    
	const fvPatchField<vector>& Up =
            patch().lookupPatchField<volVectorField, vector>(UName_);

        const basicThermo& thermo =
            db().lookupObject<basicThermo>("thermophysicalProperties");

        const volScalarField Cp = thermo.Cp();
	const volScalarField Cv = thermo.Cv();
	
	const vectorField& normalVector = patch().nf();

        const fvPatchField<scalar>& Cpp =
            patch().patchField<volScalarField, scalar>(Cp);
        
	const fvPatchField<scalar>& Cvp =
            patch().patchField<volScalarField, scalar>(Cv);	

	refValue() = fixedValue_; 
	    
	forAll(Tp, patchI)
	{
	    scalar Map = mag(normalVector[patchI] & Up[patchI])
	        /sqrt((Cpp[patchI]/Cvp[patchI])*(Cpp[patchI]-Cvp[patchI])*Tp[patchI]);
	    
	    /*
	    valueFraction()[patchI] = min( 
	      pos(Map-1.0) ? 0.0 : 1.0,
	      pos(normalVector[patchI] & Up[patchI])
	    );
	    */
	    valueFraction()[patchI] = pos(Map-1.0) ? 0.0 : 1.0;
	}

    mixedFvPatchScalarField::updateCoeffs();
}


void
mixedSubsonicSupersonicOutletFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("T") << TName_ << token::END_STATEMENT << nl;
    os.writeKeyword("U") << UName_ << token::END_STATEMENT << nl;
#if (OPENFOAM >= 1912)
    fixedValue_.writeEntry("fixedValue", os);
    this->writeEntry("value", os);
#else
    writeEntry(os, "fixedValue", fixedValue_);
    writeEntry(os, "value", *this);
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    mixedSubsonicSupersonicOutletFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
