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
    Foam::totalPressureDirectedInletVelocity

Description
    - should be used with zero gradient for pressure
    - sets the fixed value condition for the velocity, the velocity magnitude is 
      either 0 for p > p0, otherwise it is calculated using isentropic
      relations from p/p0

Example

    inlet
    {
        type            totalPressureDirectedInletVelocity;
        gamma           1.4;
        p0              uniform 100000;
        inletDirection  uniform (1 0 0);
        value           uniform 100000;
    } 

\*---------------------------------------------------------------------------*/

#ifndef totalPressureDirectedInletVelocityFvPatchVectorFields_H
#define totalPressureDirectedInletVelocityFvPatchVectorFields_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class totalPressureDirectedInletVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class totalPressureDirectedInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

    //- Name of the compressibility field used to calculate the wave speed
    word psiName_;

    //- Heat capacity ratio
    scalar gamma_;
    
    //- Total pressure
    scalarField p0_;

    //- Inlet direction
    vectorField inletDir_;

public:

    //- Runtime type information
    TypeName("totalPressureDirectedInletVelocity");


    // Constructors

        //- Construct from patch and internal field
        totalPressureDirectedInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        totalPressureDirectedInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given totalPressureDirectedInletVelocityFvPatchVectorField
        //  onto a new patch
        totalPressureDirectedInletVelocityFvPatchVectorField
        (
            const totalPressureDirectedInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        totalPressureDirectedInletVelocityFvPatchVectorField
        (
            const totalPressureDirectedInletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new totalPressureDirectedInletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        totalPressureDirectedInletVelocityFvPatchVectorField
        (
            const totalPressureDirectedInletVelocityFvPatchVectorField&,
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
                new totalPressureDirectedInletVelocityFvPatchVectorField(*this, iF)
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
