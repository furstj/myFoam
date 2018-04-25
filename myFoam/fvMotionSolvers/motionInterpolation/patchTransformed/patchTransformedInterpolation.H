/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::patchTransformedInterpolation

Description
    Interpolation of cell-based displacements to the points with additional
    correction of patch-adjacent cells.

    Snappy meshes often have moderately non-orthogonal faces next to patches.
    If these partches are moved, these already quite distorted faces quickly
    violate the quality criteria.

    This method corrects patch-adjacent cells are by calculating the average
    translation, rotation and expansion of points connected to the patch. This
    transformation is then applied to all the points in the cell that are not
    on the patch. This does a better job of preserving the cell shape on the
    patch than pure interpolation, though it can simply move the problem into
    the next layer of cells.

    The user needs only to specify the patches on which this correction is
    performed:

    \verbatim
        interpolation patchTransformed (movingPatch1 movingPatch2);
    \endverbatim

    Note that this method is not defined for scalar displacement fields.

SourceFiles
    patchTransformedInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef patchTransformedInterpolation_H
#define patchTransformedInterpolation_H

#include "motionInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class patchTransformedInterpolation Declaration
\*---------------------------------------------------------------------------*/

class patchTransformedInterpolation
:
    public motionInterpolation
{
    // Private data

        //- Patches on which to enforce transformation on adjacent cells
        const labelList patches_;


    // Private member functions

        //- Get patches from the input stream
        labelList getPatches(Istream& entry) const;


public:

    //- Runtime type information
    TypeName("patchTransformed");


    // Constructors

        //- Construct from an fvMesh and an Istream
        patchTransformedInterpolation
        (
            const fvMesh& mesh,
            Istream& entry
        );


    //- Destructor
    virtual ~patchTransformedInterpolation();


    // Member Functions

        //- Interpolate the given scalar cell displacement
        virtual void interpolate
        (
            const volScalarField&,
            pointScalarField&
        ) const;

        //- Interpolate the given vector cell displacement
        virtual void interpolate
        (
            const volVectorField&,
            pointVectorField&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
