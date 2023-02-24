/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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
    Foam::CoolPropGas

Description
    CoolPropGas wrapper for CoolProp library. 
    
    IMPORTANT NOTE:  The EOS  returns complete properties (i.e. ideal gas + departure). 
    This is different from "standard" OpenFOAM thermodynamics. Therefore it has to be
    combined with CoolProp thermo.

SourceFiles
    CoolPropGasI.H
    CoolPropGas.C

\*---------------------------------------------------------------------------*/

#ifndef CoolPropGas_H
#define CoolPropGas_H

#include "autoPtr.H"
#include "specie.H"

#include "CoolProp.h"
#include "AbstractState.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class CoolPropGas;

template<class Specie>
inline CoolPropGas<Specie> operator+
(
    const CoolPropGas<Specie>&,
    const CoolPropGas<Specie>&
);

template<class Specie>
inline CoolPropGas<Specie> operator*
(
    const scalar,
    const CoolPropGas<Specie>&
);

template<class Specie>
inline CoolPropGas<Specie> operator==
(
    const CoolPropGas<Specie>&,
    const CoolPropGas<Specie>&
);


template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const CoolPropGas<Specie>&
);



/*---------------------------------------------------------------------------*\
                           Class CoolPropGas Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class CoolPropGas
:
    public Specie
{
    // Private data
    shared_ptr<CoolProp::AbstractState> state_;

public:

    // Constructors

        //- Construct from components
        inline CoolPropGas
        (
            const Specie& sp,
            const std::string& backend, 
            const std::string& fluid_names
        );

        //- Construct from dictionary
        CoolPropGas(const dictionary& dict);

        //- Construct as named copy
        inline CoolPropGas(const word& name, const CoolPropGas&);

        //- Construct and return a clone
        inline autoPtr<CoolPropGas> clone() const;

        // Selector from dictionary
        inline static autoPtr<CoolPropGas> New
        (
            const dictionary& dict
        );


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "CoolPropGas<" + word(Specie::typeName_()) + '>';
        }

            shared_ptr<CoolProp::AbstractState> state() const
            {
                return state_;
            } 
        // Fundamental properties


            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return enthalpy [J/kg]
            inline scalar H(const scalar p, const scalar T) const;

            //- Return internal energy [J/kg]
            inline scalar E(const scalar p, const scalar T) const;

            //- Return Cp [J/(kg K)]
            inline scalar Cp(scalar p, scalar T) const;

            //- Return Cv [J/(kg K)]
            inline scalar Cv(scalar p, scalar T) const;

            //- Return entropy [J/(kg K)]
            inline scalar S(const scalar p, const scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor [-]
            inline scalar Z(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K)]
            inline scalar CpMCv(scalar p, scalar T) const;


    // IO

            //- Write to Ostream
            void write(Ostream& os) const;

    // Member operators

        inline void operator+=(const CoolPropGas&);
        inline void operator*=(const scalar);


    // Friend operators

        friend CoolPropGas operator+ <Specie>
        (
            const CoolPropGas&,
            const CoolPropGas&
        );

        friend CoolPropGas operator* <Specie>
        (
            const scalar s,
            const CoolPropGas&
        );

        friend CoolPropGas operator== <Specie>
        (
            const CoolPropGas&,
            const CoolPropGas&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const CoolPropGas&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CoolPropGasI.H"

#ifdef NoRepository
    #include "CoolPropGas.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
