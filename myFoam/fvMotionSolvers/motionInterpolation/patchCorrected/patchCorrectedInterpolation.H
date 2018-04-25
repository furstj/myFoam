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
    Foam::patchCorrectedInterpolation

Description
    Interpolation of cell-based displacements to the points with additional
    correction for interpolation inconsistency on patches.

    The default interpolation method interpolates from the cells to all points
    except those on boundaries with value boundary conditions. The discrepancy
    across the boundary cell can cause shearing and inversion if the cells are
    of very high aspect ratio.

    This method applies the default interpolation to *all* points, including
    those on value boundaries. The difference between the interpolated value on
    the boundary and the desired boundary condition is then propagated into the
    mesh with a wave. Contributions from different patches are inverse-distance
    weighted, and the result is added to the default interpolation. The result
    of this is that thin boundary cells are maintained much more accurately for
    non-uniform patch displacements.

    The user must specify the patch groups from which to propagate the motion.
    Ideally, these groups will be opposing; i.e. one group with the aerofoil,
    and one with the far field:

    \verbatim
        interpolation patchCorrected
        (
            (aerofoilUpper aerofoilLower)
            (farField)
        );
    \endverbatim

SourceFiles
    patchCorrectedInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef patchCorrectedInterpolation_H
#define patchCorrectedInterpolation_H

#include "motionInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class patchCorrectedInterpolation Declaration
\*---------------------------------------------------------------------------*/

class patchCorrectedInterpolation
:
    public motionInterpolation
{
    // Private data

        //- Patch groups from which to propagate the displacement
        const labelListList patchGroups_;


    // Private member functions

        //- Get patch groups from the input stream
        labelListList getPatchGroups(Istream& entry) const;

        //- Interpolate the given cell displacement
        template <class Type>
        void interpolateType
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            GeometricField<Type, pointPatchField, pointMesh>&
        ) const;

        //- Interpolate patch data by inverse distance weighting
        template <class Type>
        void interpolateDataFromPatchGroups
        (
            GeometricField<Type, pointPatchField, pointMesh>&
        ) const;

        //- Propagate data from a number of patches into the field
        template <class Type>
        void propagateDataFromPatchGroup
        (
            const label,
            pointScalarField&,
            GeometricField<Type, pointPatchField, pointMesh>&
        ) const;


public:

    //- Runtime type information
    TypeName("patchCorrected");


    // Constructors

        //- Construct from an fvMesh and an Istream
        patchCorrectedInterpolation
        (
            const fvMesh& mesh,
            Istream& entry
        );


    //- Destructor
    virtual ~patchCorrectedInterpolation();


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

#include "patchCorrectedInterpolationTemplates.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
