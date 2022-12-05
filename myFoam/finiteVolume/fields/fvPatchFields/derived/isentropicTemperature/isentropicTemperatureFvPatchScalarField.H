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
    Foam::isentropicTemperatureFvPatchScalarField

Group
    grpInletBoundaryConditions grpOutletBoundaryConditions

Description
    This boundary condition provides a total temperature condition calculated using real gas relations.

Usage
    \table
        Property     | Description             | Required    | Default value
        p            | Pressure field name     | no          | p
        U            | Velocity field name     | no          | U
        phi          | Flux field name         | no          | phi
        T0           | total temperature       | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            isentropicTemperature;
        gamma           1.4;
        T0              uniform 800;
    }
    \endverbatim

SourceFiles
    isentropicTemperatureFvPatchScalarField.C

See also
    Foam::fixedValueFvPatchField

\*---------------------------------------------------------------------------*/

#ifndef isentropicTemperatureFvPatchScalarField_H
#define isentropicTemperatureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class isentropicTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class isentropicTemperatureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Name of the pressure field
        word pName_;

        //- Name of the velocity field
        word UName_;

        //- Name of the flux transporting the field
        word phiName_;

        //- Total temperature
        scalarField T0_;


public:

    //- Runtime type information
    TypeName("isentropicTemperature");


    // Constructors

        //- Construct from patch and internal field
        isentropicTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        isentropicTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given isentropicTemperatureFvPatchScalarField
        //  onto a new patch
        isentropicTemperatureFvPatchScalarField
        (
            const isentropicTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        isentropicTemperatureFvPatchScalarField
        (
            const isentropicTemperatureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new isentropicTemperatureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        isentropicTemperatureFvPatchScalarField
        (
            const isentropicTemperatureFvPatchScalarField&,
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
                new isentropicTemperatureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the total pressure
            const scalarField& T0() const
            {
                return T0_;
            }

            //- Return reference to the total pressure to allow adjustment
            scalarField& T0()
            {
                return T0_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

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
