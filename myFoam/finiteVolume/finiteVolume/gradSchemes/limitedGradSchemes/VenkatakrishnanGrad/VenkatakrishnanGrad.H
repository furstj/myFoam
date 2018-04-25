/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::fv::VenkatakrishnanGrad

Description
    VenkatakrishnanGrad gradient scheme applied to a runTime selected base gradient
    scheme.

    The scalar limiter based on limiting the extrapolated face values
    between the maximum and minumum cell and cell neighbour values and is
    applied to all components of the gradient.

SourceFiles
    VenkatakrishnanGrad.C

\*---------------------------------------------------------------------------*/

#ifndef VenkatakrishnanGrad_H
#define VenkatakrishnanGrad_H

#include "gradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class VenkatakrishnanGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class VenkatakrishnanGrad
:
    public fv::gradScheme<Type>
{
    // Private Data

        tmp<fv::gradScheme<Type> > basicGradScheme_;

        //- Limiter coefficient
        const scalar k_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        VenkatakrishnanGrad(const VenkatakrishnanGrad&);

        //- Disallow default bitwise assignment
        void operator=(const VenkatakrishnanGrad&);


public:

    //- RunTime type information
    TypeName("Venkatakrishnan");


    // Constructors

        //- Construct from mesh and schemeData
        VenkatakrishnanGrad(const fvMesh& mesh, Istream& schemeData)
        :
            gradScheme<Type>(mesh),
            basicGradScheme_(fv::gradScheme<Type>::New(mesh, schemeData)),
            k_(readScalar(schemeData))
        {
            if (k_ < 0)
            {
                FatalIOErrorIn
                (
                    "VenkatakrishnanGrad(const fvMesh&, Istream& schemeData)",
                    schemeData
                )   << "coefficient = " << k_
                    << " should be >= 0"
                    << exit(FatalIOError);
            }
        }


    // Member Functions

        static inline void limitFace
        (
            Type& limiter,
	    const scalar& cellVolume,
            const Type& maxDelta,
            const Type& minDelta,
            const Type& extrapolate
        );

        //- Return the gradient of the given field to the gradScheme::grad
        //  for optional caching
        virtual tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
        > calcGrad
        (
            const GeometricField<Type, fvPatchField, volMesh>& vsf,
            const word& name
        ) const;
};


// * * * * * * * * * * * * Inline Member Function  * * * * * * * * * * * * * //

template<>
inline void VenkatakrishnanGrad<scalar>::limitFace
(
    scalar& limiter,
    const scalar& epsilonSquare,
    const scalar& maxDelta,
    const scalar& minDelta,
    const scalar& extrapolate
)
{

    if (extrapolate > 0)
    {
        limiter = min(limiter,
        (
            (sqr(maxDelta) + epsilonSquare)*extrapolate
            + 2*sqr(extrapolate)*maxDelta
        )/
        stabilise(
            extrapolate*
            (
                sqr(maxDelta)
                + 2.0*sqr(extrapolate)
                + maxDelta*extrapolate
                + epsilonSquare
            ),
            SMALL
        )
        );
    }
    else if (extrapolate < 0)
    {
        limiter = min( limiter,
        (
            (sqr(minDelta) + epsilonSquare)*extrapolate
            + 2*sqr(extrapolate)*minDelta
        )/
        stabilise
        (
            extrapolate*
            (
                sqr(minDelta)
                + 2.0*sqr(extrapolate) + minDelta*extrapolate
                + epsilonSquare
            ),
            SMALL
        )
        );
    }

}


template<class Type>
inline void VenkatakrishnanGrad<Type>::limitFace
(
    Type& limiter,
    const scalar& epsilonSquare,
    const Type& maxDelta,
    const Type& minDelta,
    const Type& extrapolate
)
{
    for (direction cmpt=0; cmpt<Type::nComponents; cmpt++)
    {
        VenkatakrishnanGrad<scalar>::limitFace
        (
            limiter.component(cmpt),
	    epsilonSquare,
            maxDelta.component(cmpt),
            minDelta.component(cmpt),
            extrapolate.component(cmpt)
        );
    }
}


// * * * * * * * * Template Member Function Specialisations  * * * * * * * * //

template<>
tmp<volVectorField> VenkatakrishnanGrad<scalar>::calcGrad
(
    const volScalarField& vsf,
    const word& name
) const;


template<>
tmp<volTensorField> VenkatakrishnanGrad<vector>::calcGrad
(
    const volVectorField& vsf,
    const word& name
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
