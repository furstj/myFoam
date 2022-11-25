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
    Foam::isentropicInletPressure

Description
    isentropicInletPressure implements a boundary condition for p with
    given total temperature.

    The static pressure is extrapolated from the internall cell, i.e. it is
    equivalent to zeroGradient condition. However, it stores the value of 
    (isentropic) total pressure for boundary conditions for U and T. 

    The boundary condition is compatible with real gas model with should be used in
    cooperation with isentropicInletVelocity and isentropicInletTemperature.

See also
    Foam::mixedFvPatchScalarField

SourceFiles
    isentropicInletPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef isentropicInletPressureFvPatchScalarField_H
#define isentropicInletPressureFvPatchScalarField_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class isentropicInletPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class isentropicInletPressureFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

    //- Total pressure
    scalarField p0_;

public:

    //- Runtime type information
    TypeName("isentropicInletPressure");


    // Constructors

        //- Construct from patch and internal field
        isentropicInletPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        isentropicInletPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  isentropicInletPressureFvPatchScalarField
        //  onto a new patch
        isentropicInletPressureFvPatchScalarField
        (
            const isentropicInletPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        isentropicInletPressureFvPatchScalarField
        (
            const isentropicInletPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new isentropicInletPressureFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        isentropicInletPressureFvPatchScalarField
        (
            const isentropicInletPressureFvPatchScalarField&,
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
                new isentropicInletPressureFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Attributes

            //- Return true: this patch field is altered by assignment
            virtual bool assignable() const
            {
                return true;
            }


        // Access

	    const scalarField& p0() const
            {
                return p0_;
            }

	    scalarField& p0() 
            {
                return p0_;
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
