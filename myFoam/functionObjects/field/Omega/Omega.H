/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::functionObjects::Omega

Group
    grpFieldFunctionObjects

Description
    Computes the vortex identification criterion Omega

    The vortex identification criterion Omega is defined as 

    \f[
        \Omega = \frac{W^2}{S^2 + W^2 + \epsilon}
    \f]

    where
    \f{align*}{
        W &= \sqrt{2 \mathbf{W} : \mathbf{W}}, \\
        S &= \sqrt{2 \mathbf{S} : \mathbf{S}}, \\
        \mathbf{W} &= \frac{1}{2} \left( \nabla \mathbf{U} - \nabla \mathbf{U}^T \right), \\
        \mathbf{S} &= \frac{1}{2} \left( \nabla \mathbf{U} + \nabla \mathbf{U}^T \right), \\
        \epsilon &= 10^{-10}
    \f}

    The vortex is identified as \f$\Omega > 0.5\f$. Typical threshold value is 0.52.

    The criterion is described in
      -  Liu, C. Q., Wang, Y. Q., Yang, Y., & Duan, Z. W. (2016). New omega vortex identification method. 
      Science China: Physics, Mechanics and Astronomy, 59(8), 1–9. https://doi.org/10.1007/S11433-016-0022-6

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    Omega
    {
        // Mandatory entries (unmodifiable)
        type        Omega;
        libs        (fieldFunctionObjects);
    }
    \endverbatim


See also
    - Foam::functionObject
    - Foam::functionObjects::fieldExpression
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::Q

SourceFiles
    Q.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_Omega_H
#define functionObjects_Omega_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class Omega Declaration
\*---------------------------------------------------------------------------*/

class Omega
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the Omega field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("Omega");


    // Constructors

        //- Construct from Time and dictionary
        Omega
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        Omega(const Omega&) = delete;

        //- No copy assignment
        void operator=(const Omega&) = delete;


    //- Destructor
    virtual ~Omega() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
