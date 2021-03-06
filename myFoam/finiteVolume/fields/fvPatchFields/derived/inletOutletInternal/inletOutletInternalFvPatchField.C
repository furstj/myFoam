/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "inletOutletInternalFvPatchField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::inletOutletInternalFvPatchField<Type>::inletOutletInternalFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    phiName_("phi")
{
    this->refValue() = pTraits<Type>::zero;
    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 0.0;
}


template<class Type>
Foam::inletOutletInternalFvPatchField<Type>::inletOutletInternalFvPatchField
(
    const inletOutletInternalFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper),
    phiName_(ptf.phiName_)
{}


template<class Type>
Foam::inletOutletInternalFvPatchField<Type>::inletOutletInternalFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi"))
{
    this->refValue() = Field<Type>("inletValue", dict, p.size());

    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<Type>::operator=(this->refValue());
    }

    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 0.0;
}


template<class Type>
Foam::inletOutletInternalFvPatchField<Type>::inletOutletInternalFvPatchField
(
    const inletOutletInternalFvPatchField<Type>& ptf
)
:
    mixedFvPatchField<Type>(ptf),
    phiName_(ptf.phiName_)
{}


template<class Type>
Foam::inletOutletInternalFvPatchField<Type>::inletOutletInternalFvPatchField
(
    const inletOutletInternalFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF),
    phiName_(ptf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::inletOutletInternalFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    /*const Field<scalar>& phip =
        this->patch().template lookupPatchField<surfaceScalarField, scalar>
        (
            phiName_
        );

    this->valueFraction() = 1.0 - pos(phip);
    */

    const vectorField& normal = this->patch().nf();
    
    const fvPatchField<vector>& Up =
      this->patch().template lookupPatchField<volVectorField, vector>("U");

    const vectorField& Uint = Up.internalField();

    this->valueFraction() = 1.0 - pos( normal & Uint );

    mixedFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::inletOutletInternalFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    if (phiName_ != "phi")
    {
        os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    }
#if (OPENFOAM >= 1812)
    const Field<Type>& refValue_ = this->refValue();
    refValue_.writeEntry("inletValue", os);
    this->writeEntry("value", os);
#else
    writeEntry(os, "inletValue", this->refValue());
    writeEntry(os, "value", *this);
#endif
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::inletOutletInternalFvPatchField<Type>::operator=
(
    const fvPatchField<Type>& ptf
)
{
    fvPatchField<Type>::operator=
    (
        this->valueFraction()*this->refValue()
        + (1 - this->valueFraction())*ptf
    );
}


// ************************************************************************* //
