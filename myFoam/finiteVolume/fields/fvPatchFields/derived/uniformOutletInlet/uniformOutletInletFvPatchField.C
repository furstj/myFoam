/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
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

#include "uniformOutletInletFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::uniformOutletInletFvPatchField<Type>::uniformOutletInletFvPatchField
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
Foam::uniformOutletInletFvPatchField<Type>::uniformOutletInletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    uniformOutletValue_(Function1<Type>::New("uniformOutletValue", dict))
{
    this->refValue() =
        uniformOutletValue_->value(this->db().time().timeOutputValue());

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
Foam::uniformOutletInletFvPatchField<Type>::uniformOutletInletFvPatchField
(
    const uniformOutletInletFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(p, iF),  // Don't map
    phiName_(ptf.phiName_),
#if OPENFOAM_PLUS >= 1812
    uniformOutletValue_(ptf.uniformOutletValue_.clone())
#else
    uniformOutletValue_(ptf.uniformOutletValue_, false)
#endif
{
    this->patchType() = ptf.patchType();

    // Evaluate refValue since not mapped
    this->refValue() =
        uniformOutletValue_->value(this->db().time().timeOutputValue());

    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 0.0;

    // Initialize the patch value to the refValue
    fvPatchField<Type>::operator=(this->refValue());

    this->map(ptf, mapper);
}


template<class Type>
Foam::uniformOutletInletFvPatchField<Type>::uniformOutletInletFvPatchField
(
    const uniformOutletInletFvPatchField<Type>& ptf
)
:
    mixedFvPatchField<Type>(ptf),
    phiName_(ptf.phiName_),
#if OPENFOAM_PLUS >= 1812
    uniformOutletValue_(ptf.uniformOutletValue_.clone())
#else
    uniformOutletValue_(ptf.uniformOutletValue_, false)
#endif	
{}


template<class Type>
Foam::uniformOutletInletFvPatchField<Type>::uniformOutletInletFvPatchField
(
    const uniformOutletInletFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF),
    phiName_(ptf.phiName_),
#if OPENFOAM_PLUS >= 1812
    uniformOutletValue_(ptf.uniformOutletValue_.clone())
#else
    uniformOutletValue_(ptf.uniformOutletValue_, false)
#endif
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::uniformOutletInletFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    this->refValue() =
        uniformOutletValue_->value(this->db().time().timeOutputValue());

    const Field<scalar>& phip =
        this->patch().template lookupPatchField<surfaceScalarField, scalar>
        (
            phiName_
        );

    this->valueFraction() = pos(phip);

    mixedFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::uniformOutletInletFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    if (phiName_ != "phi")
    {
        os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    }
    this->uniformOutletValue_->writeData(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::uniformOutletInletFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchField<Type>::autoMap(m);

    // Override
    this->refValue() =
        uniformOutletValue_->value(this->db().time().timeOutputValue());
}


template<class Type>
void Foam::uniformOutletInletFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<Type>::rmap(ptf, addr);

    // Override
    this->refValue() =
        uniformOutletValue_->value(this->db().time().timeOutputValue());
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::uniformOutletInletFvPatchField<Type>::operator=
(
    const fvPatchField<Type>& ptf
)
{
    fvPatchField<Type>::operator=
    (
     (1 - this->valueFraction())*this->refValue()
        + this->valueFraction()*ptf
    );
}


// ************************************************************************* //
