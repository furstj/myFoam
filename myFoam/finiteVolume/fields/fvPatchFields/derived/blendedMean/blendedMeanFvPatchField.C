/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2019 OpenFOAM Foundation
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

#include "blendedMeanFvPatchField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::blendedMeanFvPatchField<Type>::blendedMeanFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    meanValue_(),
    blendFactor_()
{}


template<class Type>
Foam::blendedMeanFvPatchField<Type>::blendedMeanFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    meanValue_(Function1<Type>::New("meanValue", dict)),
    blendFactor_(Function1<scalar>::New("blendFactor", dict))
{}


template<class Type>
Foam::blendedMeanFvPatchField<Type>::blendedMeanFvPatchField
(
    const blendedMeanFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
#if (OPENFOAM >= 1812)
    meanValue_(ptf.meanValue_.clone()),
    blendFactor_(ptf.blendFactor_.clone())
#else
    meanValue_(ptf.meanValue_, false),
    blendFactor_(ptf.blendFactor_, false)
#endif
{}


template<class Type>
Foam::blendedMeanFvPatchField<Type>::blendedMeanFvPatchField
(
    const blendedMeanFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
#if (OPENFOAM >= 1812)
    meanValue_(ptf.meanValue_.clone()),
    blendFactor_(ptf.blendFactor_.clone())
#else
    meanValue_(ptf.meanValue_, false),
    blendFactor_(ptf.blendFactor_, false)
#endif
{}


template<class Type>
Foam::blendedMeanFvPatchField<Type>::blendedMeanFvPatchField
(
    const blendedMeanFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
#if (OPENFOAM >= 1812)
    meanValue_(ptf.meanValue_.clone()),
    blendFactor_(ptf.blendFactor_.clone())
#else
    meanValue_(ptf.meanValue_, false),
    blendFactor_(ptf.blendFactor_, false)
#endif
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::blendedMeanFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();
    Type meanValue = meanValue_->value(t);
    scalar blendFactor = blendFactor_->value(t);

    Field<Type> newValues(this->patchInternalField());

    Type meanValuePsi =
        gSum(this->patch().magSf()*newValues)
       /gSum(this->patch().magSf());

    if (mag(meanValue) > SMALL && mag(meanValuePsi)/mag(meanValue) > 0.5)
    {
        newValues *= mag(meanValue)/mag(meanValuePsi);
    }
    else
    {
        newValues += (meanValue - meanValuePsi);
    }

    this->operator==(blendFactor*meanValue + (1-blendFactor)*newValues);

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::blendedMeanFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    #if (OPENFOAM >= 1812)
    meanValue_->writeData(os);
    blendFactor_->writeData(os);
    fvPatchField<Type>::writeValueEntry(os);
    #else 
    writeEntry(os, meanValue_());
    writeEntry(os, blendFactor_());
    writeEntry(os, "value", *this);
    #endif
}


// ************************************************************************* //
