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
    Foam::subsonicInletTotal

Description
    subsonicInletTotal implements a boundary condition for U with
    given total temperature and given direction.  

    The boundary value of velosity is computed using Riemann
    invariant going along the outward directed characteristics 
    (see e.g. Blazek)

    The boundary condition is valid only for compressible ideal gas and it is
    compatible with totalTemperature for T and totalPressure for p.

See also
    Foam::mixedFvPatchVectorField

SourceFiles
    subsonicInletTotalFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef subsonicInletTotalFvPatchVectorField_H
#define subsonicInletTotalFvPatchVectorField_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class subsonicInletTotalFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class subsonicInletTotalFvPatchVectorField
:
    public mixedFvPatchVectorField
{
    // Private data

        //- Temperature field name
        word TName_;

       
        //- Inlet direction
        vectorField inletDir_;


public:

    //- Runtime type information
    TypeName("subsonicInletTotal");


    // Constructors

        //- Construct from patch and internal field
        subsonicInletTotalFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        subsonicInletTotalFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  subsonicInletTotalFvPatchVectorField
        //  onto a new patch
        subsonicInletTotalFvPatchVectorField
        (
            const subsonicInletTotalFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        subsonicInletTotalFvPatchVectorField
        (
            const subsonicInletTotalFvPatchVectorField&
        ) = delete;

        //- Construct as copy setting internal field reference
        subsonicInletTotalFvPatchVectorField
        (
            const subsonicInletTotalFvPatchVectorField&,
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
                new subsonicInletTotalFvPatchVectorField
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

        virtual void operator=(const fvPatchField<vector>& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
