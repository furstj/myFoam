/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::CoolPropTransport

Description
    Transport package using CoolProp library.

    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    CoolPropTransportI.H
    CoolPropTransport.C

\*---------------------------------------------------------------------------*/

#ifndef CoolPropTransport_H
#define CoolPropTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class CoolPropTransport;

template<class Thermo>
inline CoolPropTransport<Thermo> operator+
(
    const CoolPropTransport<Thermo>&,
    const CoolPropTransport<Thermo>&
);

template<class Thermo>
inline CoolPropTransport<Thermo> operator*
(
    const scalar,
    const CoolPropTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const CoolPropTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                    Class CoolPropTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class CoolPropTransport
:
    public Thermo
{
    // Private data

public:

    // Constructors

        //- Construct from components
        inline CoolPropTransport
        (
            const Thermo& t
        );

        //- Construct as named copy
        inline CoolPropTransport(const word&, const CoolPropTransport&);

        //- Construct from dictionary
        CoolPropTransport(const dictionary& dict);

        //- Construct from base thermo and dictionary
        CoolPropTransport(const Thermo& t, const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<CoolPropTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<CoolPropTransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "CoolProp<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        // inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline void operator=(const CoolPropTransport&);

        inline void operator+=(const CoolPropTransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend CoolPropTransport operator+ <Thermo>
        (
            const CoolPropTransport&,
            const CoolPropTransport&
        );

        friend CoolPropTransport operator* <Thermo>
        (
            const scalar,
            const CoolPropTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const CoolPropTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CoolPropTransportI.H"

#ifdef NoRepository
    #include "CoolPropTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
