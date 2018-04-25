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

Class
    Foam::mixedSubsonicSupersonicOutletFvPatchScalarField

Description
    Foam::mixedSubsonicSupersonicOutletFvPatchScalarField
    
    BC switches between fixedValue and zeroGradient depending on the
    Mach number.
    If the outflow Mach number is subsonic, the boundary condition
    uses the specified value. Otherwhise it is zeroGrandient.
    
    Setting up the BC:
        type            mixedSubsonicSupersonicOutlet;
        fixedValue      uniform 101148;
	value		101148;    
    
Author
    Sebastian Saegeler  All rights reserved.

SourceFiles
    mixedSubsonicSupersonicOutletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef mixedSubsonicSupersonicOutletFvPatchScalarField_H
#define mixedSubsonicSupersonicOutletFvPatchScalarField_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class mixedSubsonicSupersonicOutletFvPatch Declaration
\*---------------------------------------------------------------------------*/

class mixedSubsonicSupersonicOutletFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Name of static Temperature field
        word TName_;
	
	//- fixed value if outflow is subsonic
	scalarField fixedValue_;

        //- Name of Velocity field
        word UName_;
	

public:

    //- Runtime type information
    TypeName("mixedSubsonicSupersonicOutlet");


    // Constructors

        //- Construct from patch and internal field
        mixedSubsonicSupersonicOutletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mixedSubsonicSupersonicOutletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  mixedSubsonicSupersonicOutletFvPatchScalarField
        //  onto a new patch
        mixedSubsonicSupersonicOutletFvPatchScalarField
        (
            const mixedSubsonicSupersonicOutletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mixedSubsonicSupersonicOutletFvPatchScalarField
        (
            const mixedSubsonicSupersonicOutletFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new mixedSubsonicSupersonicOutletFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        mixedSubsonicSupersonicOutletFvPatchScalarField
        (
            const mixedSubsonicSupersonicOutletFvPatchScalarField&,
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
                new mixedSubsonicSupersonicOutletFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Access

           //- Return the fixed value
            const scalarField& fixedValue() const
            {
                return fixedValue_;
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
