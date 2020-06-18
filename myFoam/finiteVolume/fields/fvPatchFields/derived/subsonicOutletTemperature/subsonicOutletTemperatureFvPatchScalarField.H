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
    Foam::subsonicOutletVelocity

Description
    subsonicOutletTemperature implements an outlet boundary condition for T based 
    characteristics, see e.g. Blazek.

    Let p denotes the outlet value and i is the internal field value. Then
    the boundary condition for u reads:

    rhop = rhoi + (pi - pp)/ci^2
    
    and then (using  ideal gas EOS)
 
    Tp = pp/(rhop*R)

    The boundary condition is valid only at subsonic inlet with prescribed pressure.

\*---------------------------------------------------------------------------*/

#ifndef subsonicOutletTemperatureFvPatchScalarFields_H
#define subsonicOutletTemperatureFvPatchScalarFields_H

#include "fvPatchField.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class subsonicOutletVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class subsonicOutletTemperatureFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

    //- Variable storing pressure
    word pName_;

    //- Variable storing density
    word rhoName_;

    //- Variable storing compressibility
    word psiName_;

public:

    //- Runtime type information
    TypeName("subsonicOutletTemperature");


    // Constructors

        //- Construct from patch and internal field
        subsonicOutletTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        subsonicOutletTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given subsonicOutletTemperatureFvPatchScalarField
        //  onto a new patch
        subsonicOutletTemperatureFvPatchScalarField
        (
            const subsonicOutletTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        subsonicOutletTemperatureFvPatchScalarField
        (
            const subsonicOutletTemperatureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new subsonicOutletTemperatureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        subsonicOutletTemperatureFvPatchScalarField
        (
            const subsonicOutletTemperatureFvPatchScalarField&,
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
                new subsonicOutletTemperatureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

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

            //- Return the name of rho
            const word& rhoName() const
            {
                return rhoName_;
            }

            //- Return reference to the name of rho to allow adjustment
            word& rhoName()
            {
                return rhoName_;
            }

            //- Return the name of psi
            const word& psiName() const
            {
                return psiName_;
            }

            //- Return reference to the name of T to allow adjustment
            word& psiName()
            {
                return psiName_;
            }

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
