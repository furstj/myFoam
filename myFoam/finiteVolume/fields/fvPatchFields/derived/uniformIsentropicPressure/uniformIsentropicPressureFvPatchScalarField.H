/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::uniformIsentropicPressureFvPatchScalarField

Group
    grpInletBoundaryConditions grpOutletBoundaryConditions

Description
    This boundary condition provides a total pressure condition for real gas.

Usage
    \table
        Property     | Description                | Required | Default value
        U            | Velocity field name        | no       | U
        T            | Temperature fielad name    | no       | T
        phi          | Flux field name            | no       | phi
        p0           | Total pressure as a function of time  | yes   |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            uniformIsentropicPressure;
        p0             uniform 1e5;
    }
    \endverbatim

    The \c p0 entry is specified as a Function1 type, able to describe
    time varying functions.

See also
    Foam::DataEntry
    Foam::fixedValueFvPatchField
    Foam::isentropicPressureFvPatchField

SourceFiles
    uniformIsentropicPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef uniformIsentropicPressureFvPatchScalarField_H
#define uniformIsentropicPressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class uniformIsentropicPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class uniformIsentropicPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Name of the velocity field
        word UName_;

        //- Name of the temperature field
        word TName_;

        //- Name of the flux transporting the field
        word phiName_;

        //- Total pressure
        autoPtr<Function1<scalar> > p0_;


public:

    //- Runtime type information
    TypeName("uniformIsentropicPressure");


    // Constructors

        //- Construct from patch and internal field
        uniformIsentropicPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        uniformIsentropicPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given uniformIsentropicPressureFvPatchScalarField
        //  onto a new patch
        uniformIsentropicPressureFvPatchScalarField
        (
            const uniformIsentropicPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        uniformIsentropicPressureFvPatchScalarField
        (
            const uniformIsentropicPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new uniformIsentropicPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        uniformIsentropicPressureFvPatchScalarField
        (
            const uniformIsentropicPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new uniformIsentropicPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the name of the velocity field
            const word& UName() const
            {
                return UName_;
            }

            //- Return reference to the name of the velocity field
            //  to allow adjustment
            word& UName()
            {
                return UName_;
            }

            //- Return the name of the temperature field
            const word& TName() const
            {
                return TName_;
            }

            //- Return reference to the name of the temperature field
            //  to allow adjustment
            word& TName()
            {
                return TName_;
            }

            //- Return the name of the flux field
            const word& phiName() const
            {
                return phiName_;
            }

            //- Return reference to the name of the flux field
            //  to allow adjustment
            word& phiName()
            {
                return phiName_;
            }


        // Evaluation functions

            //- Inherit updateCoeffs from fixedValueFvPatchScalarField
            using fixedValueFvPatchScalarField::updateCoeffs;

            //- Update the coefficients associated with the patch field
            //  using the given patch velocity field
            virtual void updateCoeffs(const vectorField& Up);

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
