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
    Foam::isentropicInletVelocity

Description
    isentropicInletVelocity implements a boundary condition for U with
    given total temperature and given direction or tangential velocity.  

    The boundary condition is valid for compressible real gas and it is
    compatible with isentropicTemperature for T and isentropicPressure for p.
    
    The inlet direction is either:
    - given by an optional inletDirection field, or
    - has given tangential component given by an optional tangentialVelocity, or
    - parallel to inward normal if no 'inletDirection' or 'tangentialVelocity' are given

See also
    Foam::mixedFvPatchVectorField

SourceFiles
    isentropicInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef isentropicInletVelocityFvPatchVectorField_H
#define isentropicInletVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class isentropicInletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class isentropicInletVelocityFvPatchVectorField
:
    public mixedFvPatchVectorField
{
    // Private data

        //- Pressure field name
        word pName_;

        //- Temperature field name
        word TName_;

       
        //- Inlet direction
        vectorField inletDir_;

        //- Optional tangential velocity
        //- Note: either inletDir, or tangentialVelocity can be specified, but no both!
        vectorField tangentialVelocity_;

public:

    //- Runtime type information
    TypeName("isentropicInletVelocity");


    // Constructors

        //- Construct from patch and internal field
        isentropicInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        isentropicInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  isentropicInletVelocityFvPatchVectorField
        //  onto a new patch
        isentropicInletVelocityFvPatchVectorField
        (
            const isentropicInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        isentropicInletVelocityFvPatchVectorField
        (
            const isentropicInletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new isentropicInletVelocityFvPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        isentropicInletVelocityFvPatchVectorField
        (
            const isentropicInletVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new isentropicInletVelocityFvPatchVectorField
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

            //- Return the name of p
            const word& pName() const
            {
                return pName_;
            }

            //- Return reference to the name of p to allow adjustment
            word& pName()
            {
                return pName_;
            }

            //- Return the name of T
            const word& TName() const
            {
                return TName_;
            }

            //- Return reference to the name of T to allow adjustment
            word& TName()
            {
                return TName_;
            }

            const vectorField& inletDirection() const
            {
		return inletDir_;
            }

	    vectorField& inletDirection() 
            {
                return inletDir_;
            }

    	    const vectorField& tangentialVelocity() 
            {
                return tangentialVelocity_;
            }

            //- Reset the tangential velocity
            void setTangentialVelocity(const vectorField& tangentialVelocity);

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        //virtual void operator=(const fvPatchField<vector>& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
