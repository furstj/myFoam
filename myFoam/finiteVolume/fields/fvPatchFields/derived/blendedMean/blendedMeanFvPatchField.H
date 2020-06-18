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

Class
    Foam::blendedMeanFvPatchField

Description
    This boundary condition extrapolates field to the patch using the near-cell
    values and adjusts the distribution to match the specified, optionally
    time-varying, mean value.

    The boundary value is then blended with (time-varying) uniform profile
    as

    val_f = blendFactor*meanValue + 
            (1-blendFactor)*(val_i + (meanValue-avg(val_i))

    i.e. the condition is equivalent to fixedValue if blendFactor=1 and 
    to fixedMean when blendFactor=0.

Usage
    \table
        Property     | Description             | Required    | Default value
        meanValue    | mean value Function1    | yes         |
        blendFactor  | blend factor Function1  | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fixedMean;
        meanValue       1.0;
	blendFactor     0.05;
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField
    Foam::Function1Types
    fixedMeanFvPatchField.C

SourceFiles
    blendedMeanFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef blendedMeanFvPatchField_H
#define blendedMeanFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class blendedMeanFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class blendedMeanFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private Data

        //- MeanValue value the field is adjusted to maintain
        autoPtr<Function1<Type>> meanValue_;

        //- BlendFactor blending between fixedVale and fixedmean
        autoPtr<Function1<scalar>> blendFactor_;


public:

    //- Runtime type information
    TypeName("blendedMean");


    // Constructors

        //- Construct from patch and internal field
        blendedMeanFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        blendedMeanFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given blendedMeanFvPatchField
        //  onto a new patch
        blendedMeanFvPatchField
        (
            const blendedMeanFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        blendedMeanFvPatchField
        (
            const blendedMeanFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new blendedMeanFvPatchField<Type>(*this)
            );
        }

        //- Copy constructor setting internal field reference
        blendedMeanFvPatchField
        (
            const blendedMeanFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new blendedMeanFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "blendedMeanFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
