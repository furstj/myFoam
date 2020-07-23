
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

#include "mixingPlaneFvPatchField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::mixingPlaneFvPatchField<Type>::mixingPlaneFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    origin_(Zero),
    axis_(Zero),
    order_(Zero),
    source_()
{}


template<class Type>
Foam::mixingPlaneFvPatchField<Type>::mixingPlaneFvPatchField
(
    const mixingPlaneFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    order_(ptf.order_),
    source_(ptf.source_)
{}


template<class Type>
Foam::mixingPlaneFvPatchField<Type>::mixingPlaneFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    origin_(dict.lookupOrDefault<vector>("origin", vector(0,0,0))),
    axis_(dict.lookupOrDefault<vector>("axis", vector(1,0,0))),
    order_(dict.lookupOrDefault<label>("order", 0)),
    source_(dict.lookup("source"))
{
    axis_ /= mag(axis_);
}


template<class Type>
Foam::mixingPlaneFvPatchField<Type>::mixingPlaneFvPatchField
(
    const mixingPlaneFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    order_(ptf.order_),
    source_(ptf.source_)
{}


template<class Type>
Foam::mixingPlaneFvPatchField<Type>::mixingPlaneFvPatchField
(
    const mixingPlaneFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    order_(ptf.order_),
    source_(ptf.source_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::mixingPlaneFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const GeometricField<Type, fvPatchField, volMesh>& f
        (
            dynamic_cast<const GeometricField<Type, fvPatchField, volMesh>&>
            (
                this->internalField()
            )
        );
    
    const fvPatch& p = this->patch();
    label sourcePatchID =
        p.patch().boundaryMesh().findPatchID(this->source_);
    
    if (sourcePatchID < 0)
    {
        FatalErrorInFunction
            << "Unable to find source patch " << source_
                << abort(FatalError);
    }

    const fvPatch& sourcePatch = p.boundaryMesh()[sourcePatchID];

    const fvPatchField<Type>& sourcePatchField =
        f.boundaryField()[sourcePatchID];

    // Build a system of normal equations for polynomial approximation
    label n = this->order_ + 1;
    scalarSymmetricSquareMatrix ATA(n, scalar(0));
    List<Type> ATb(n, Type(Zero));
    
    forAll(sourcePatchField, i)
    {
        vector r = sourcePatch.Cf()[i] - this->origin_;
        r -= (this->axis_ & r)*this->axis_;
        scalar w = sourcePatch.magSf()[i];
        scalar xi = mag(r);
        Type   yi = toXRTheta(sourcePatchField[i], r/mag(r));
        for (label j=0; j<n; j++)
        {
            for (label k=0; k<=j; k++)
            {
                ATA(j,k) += sqr(w)*pow(xi,j+k);
            }
            ATb[j] += sqr(w)*pow(xi,j)*yi;
        }        
    }
#if (OPENFOAM >= 1912)
    for (label i=0; i<n; i++)
        for (label j=0; j<n; j++)
            reduce(ATA(i,j), sumOp<scalar>());
#else
    reduce(ATA, sumOp<scalarSymmetricSquareMatrix>());
#endif
    reduce(ATb, sumOp<List<Type>>());

    LUsolve(ATA, ATb);
    
    Field<Type>& patchField = *this;
    forAll(patchField, i)
    {
        vector r = p.Cf()[i] - this->origin_;
        r -= (this->axis_ & r)*this->axis_;
        scalar x = mag(r);
        Type val(Zero);
        for (label j=0; j<n; j++)
        {
            val += ATb[j]*pow(x,j);
        }
        patchField[i] = fromXRTheta(val, r/mag(r));
    }
    
    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::mixingPlaneFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("origin") << origin_ << token::END_STATEMENT << nl;
    os.writeKeyword("axis") << axis_ << token::END_STATEMENT << nl;
    os.writeKeyword("order") << order_ << token::END_STATEMENT << nl;
    os.writeKeyword("source") << source_ << token::END_STATEMENT << nl;
#if (OPENFOAM >= 1912)
    this->writeEntry("value", os);
#else
    writeEntry(os, "value", *this);
#endif
}


// Specializations of to/from XRTheta

namespace Foam {

template<class Type>
Type mixingPlaneFvPatchField<Type>::toXRTheta(const Type& data, const vector& rDir) const
{
    NotImplemented;
    return data;
}

template<class Type>
Type mixingPlaneFvPatchField<Type>::fromXRTheta(const Type& data, const vector& rDir) const
{
    NotImplemented;
    return data;
}


template<>
scalar mixingPlaneFvPatchField<scalar>::toXRTheta(const scalar& data, const vector& rDir) const
{
    return data;
}

template<>
Foam::scalar mixingPlaneFvPatchField<scalar>::fromXRTheta(const scalar& data, const vector& rDir) const
{
    return data;
}


template<>
vector mixingPlaneFvPatchField<vector>::toXRTheta(const vector& data, const vector& rDir) const
{
    vector theta = axis_ ^ rDir;
    return vector(data & axis_, data & rDir, data & theta);
}

template<>
vector mixingPlaneFvPatchField<vector>::fromXRTheta(const vector& data, const vector& rDir) const
{
    vector theta = axis_ ^ rDir;
    return data[0]*axis_ + data[1]*rDir + data[2]*theta;
}

}

// ************************************************************************* //
