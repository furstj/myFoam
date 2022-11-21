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
    Foam::AungierRedlichKwongGas

Description
    AungierRedlichKwongGas gas equation of state.

SourceFiles
    AungierRedlichKwongGasI.H
    AungierRedlichKwongGas.C

\*---------------------------------------------------------------------------*/

#ifndef AungierRedlichKwongGas_H
#define AungierRedlichKwongGas_H

#include "autoPtr.H"
#include "specie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class AungierRedlichKwongGas;

template<class Specie>
inline AungierRedlichKwongGas<Specie> operator+
(
    const AungierRedlichKwongGas<Specie>&,
    const AungierRedlichKwongGas<Specie>&
);

template<class Specie>
inline AungierRedlichKwongGas<Specie> operator*
(
    const scalar,
    const AungierRedlichKwongGas<Specie>&
);

template<class Specie>
inline AungierRedlichKwongGas<Specie> operator==
(
    const AungierRedlichKwongGas<Specie>&,
    const AungierRedlichKwongGas<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const AungierRedlichKwongGas<Specie>&
);



/*---------------------------------------------------------------------------*\
                           Class AungierRedlichKwongGas Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class AungierRedlichKwongGas
:
    public Specie
{
    // Private data

        //- Critical Temperature [K]
        scalar Tc_;

        //- Critical volume [m^3/kg]
        scalar Vc_;

        //- Critical Pressure [Pa]
        scalar Pc_;

        //- Acentric factor [-]
        scalar omega_;

        //- parameter a [m^5/(kg*s^2)]]
        scalar a_;

        //- parameter b [m^3/kg]
        scalar b_;

        //- parameter b [m^3/kg]
        scalar c_;

        //- parameter n [-]
        scalar n_;

    // Private methods
    inline Foam::scalar dpdV(scalar V, scalar T, scalar alpha) const;

    inline Foam::scalar dpdT(scalar V, scalar T, scalar alpha) const;


public:

    // Constructors

        //- Construct from components
        inline AungierRedlichKwongGas
        (
            const Specie& sp,
            const scalar& Tc,
            const scalar& Vc,
            const scalar& Pc,
            const scalar& omega
        );

        //- Construct from dictionary
        AungierRedlichKwongGas(const dictionary& dict);

        //- Construct as named copy
        inline AungierRedlichKwongGas(const word& name, const AungierRedlichKwongGas&);

        //- Construct and return a clone
        inline autoPtr<AungierRedlichKwongGas> clone() const;

        // Selector from dictionary
        inline static autoPtr<AungierRedlichKwongGas> New
        (
            const dictionary& dict
        );


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "AungierRedlichKwongGas<" + word(Specie::typeName_()) + '>';
        }

        // Fundamental properties


            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return enthalpy departure [J/kg]
            inline scalar H(const scalar p, const scalar T) const;

            //- Return Cp departure [J/(kg K)]
            inline scalar Cp(scalar p, scalar T) const;

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

        inline void operator+=(const AungierRedlichKwongGas&);
        inline void operator*=(const scalar);


    // Friend operators

        friend AungierRedlichKwongGas operator+ <Specie>
        (
            const AungierRedlichKwongGas&,
            const AungierRedlichKwongGas&
        );

        friend AungierRedlichKwongGas operator* <Specie>
        (
            const scalar s,
            const AungierRedlichKwongGas&
        );

        friend AungierRedlichKwongGas operator== <Specie>
        (
            const AungierRedlichKwongGas&,
            const AungierRedlichKwongGas&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const AungierRedlichKwongGas&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "AungierRedlichKwongGasI.H"

#ifdef NoRepository
    #include "AungierRedlichKwongGas.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
