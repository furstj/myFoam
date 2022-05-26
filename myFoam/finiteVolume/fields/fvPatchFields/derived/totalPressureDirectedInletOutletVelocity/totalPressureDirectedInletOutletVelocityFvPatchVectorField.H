/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
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
    Foam::totalPressureDirectedInletOutletVelocity

Description
    - should be used with zero gradient for pressure

\*---------------------------------------------------------------------------*/

#ifndef totalPressureDirectedInletOutletVelocityFvPatchVectorFields_H
#define totalPressureDirectedInletOutletVelocityFvPatchVectorFields_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class totalPressureDirectedInletOutletVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class totalPressureDirectedInletOutletVelocityFvPatchVectorField
:
    public mixedFvPatchVectorField
{
    // Private data

    //- Name of the compressibility field used to calculate the wave speed
    word psiName_;

    //- Name of the mass flux for switching between inlet and outlet
    word phiName_;
    
    //- Heat capacity ratio
    scalar gamma_;
    
    //- Total pressure
    scalarField p0_;

    //- Inlet direction
    vectorField inletDir_;

public:

    //- Runtime type information
    TypeName("totalPressureDirectedInletOutletVelocity");


    // Constructors

        //- Construct from patch and internal field
        totalPressureDirectedInletOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        totalPressureDirectedInletOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given totalPressureDirectedInletOutletVelocityFvPatchVectorField
        //  onto a new patch
        totalPressureDirectedInletOutletVelocityFvPatchVectorField
        (
            const totalPressureDirectedInletOutletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        totalPressureDirectedInletOutletVelocityFvPatchVectorField
        (
            const totalPressureDirectedInletOutletVelocityFvPatchVectorField&
        ) = delete;
    

        //- Construct as copy setting internal field reference
        totalPressureDirectedInletOutletVelocityFvPatchVectorField
        (
            const totalPressureDirectedInletOutletVelocityFvPatchVectorField&,
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
                new totalPressureDirectedInletOutletVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            virtual const vectorField& inletDir() 
            {
                return inletDir_;
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
                const fvPatchVectorField&,
                const labelList&
            );


    virtual void updateCoeffs();
    
    
    //- Write
    virtual void write(Ostream&) const;
    
    virtual void operator=(const fvPatchField<vector>& pvf);
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
