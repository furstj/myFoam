/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::uniformOutletInletFvPatchField

Group
    grpOutletBoundaryConditions

Description
    Variant of outletInlet boundary condition with uniform outletValue.

    \heading Patch usage

    \table
        Property             | Description             | Required    | Default value
        phi                  | flux field name         | no          | phi
        uniformOutletValue   | outlet value for outlet flow | yes    |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type                uniformOutletInlet;
        phi                 phi;
        uniformOutletValue  0;
        value               uniform 0;
    }
    \endverbatim

    The mode of operation is determined by the sign of the flux across the
    patch faces.

Note
    Sign conventions:
    - negative flux (into domain): apply zero-gradient condition
    - positive flux (out of domain): apply the user-specified fixed value

SeeAlso
    Foam::inletOutletFvPatchField

SourceFiles
    uniformOutletInletFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef uniformOutletInletFvPatchField_H
#define uniformOutletInletFvPatchField_H

#include "mixedFvPatchField.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class uniformOutletInletFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class uniformOutletInletFvPatchField
:
    public mixedFvPatchField<Type>
{

protected:

    // Protected data

        //- Name of flux field
        word phiName_;

        //- Value
        autoPtr<Function1<Type> > uniformOutletValue_;


public:

    //- Runtime type information
    TypeName("uniformOutletInlet");


    // Constructors

        //- Construct from patch and internal field
        uniformOutletInletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        uniformOutletInletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given uniformOutletInletFvPatchField
        //  onto a new patch
        uniformOutletInletFvPatchField
        (
            const uniformOutletInletFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        uniformOutletInletFvPatchField
        (
            const uniformOutletInletFvPatchField<Type>&
        ) = delete;

        //- Construct as copy setting internal field reference
        uniformOutletInletFvPatchField
        (
            const uniformOutletInletFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new uniformOutletInletFvPatchField<Type>(*this, iF)
            );
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
            const fvPatchField<Type>&,
            const labelList&
        );


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const fvPatchField<Type>& pvf);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "uniformOutletInletFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
