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
    Foam::CoolPropmetaThermo

Description
    CoolProp tables based thermodynamics package templated
    into the equation of state.

    IMPORTANT NOTE: the package does not split quantities to ideal gas + departure!
    It has to be combined with CoolProp EOS which returns itself the sum
    of real gas and departure parts.

SourceFiles
    CoolPropThermoI.H
    CoolPropThermo.C

\*---------------------------------------------------------------------------*/

#ifndef CoolPropThermo_H
#define CoolPropThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class CoolPropThermo;

template<class EquationOfState>
inline CoolPropThermo<EquationOfState> operator+
(
    const CoolPropThermo<EquationOfState>&,
    const CoolPropThermo<EquationOfState>&
);

template<class EquationOfState>
inline CoolPropThermo<EquationOfState> operator*
(
    const scalar,
    const CoolPropThermo<EquationOfState>&
);

template<class EquationOfState>
inline CoolPropThermo<EquationOfState> operator==
(
    const CoolPropThermo<EquationOfState>&,
    const CoolPropThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const CoolPropThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                           Class CoolPropThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class CoolPropThermo
:
    public EquationOfState
{

public:

    // Public data
    
private:
    
    // Private data
    
    
public:
    
    // Constructors
    
    //- Construct from components
    inline CoolPropThermo
    (
        const EquationOfState& st
    );
    
    //- Construct from dictionary
    CoolPropThermo(const dictionary& dict);
    
    //- Construct as a named copy
    inline CoolPropThermo(const word&, const CoolPropThermo&);
    
    
    // Member Functions
    
    //- Return the instantiated type name
    static word typeName()
        {
            return "CoolProp<" + EquationOfState::typeName() + '>';
        }
    
    //- Limit the temperature to be in the range Tlow_ to Thigh_
    inline scalar limit(const scalar T) const;
    

    // Fundamental properties
    
    //- Heat capacity at constant pressure [J/(kg K)]
    inline scalar Cp(const scalar p, const scalar T) const;

    //- Heat capacity at constant volume [J/(kg K)]
    inline scalar Cv(const scalar p, const scalar T) const;

    //- Absolute Enthalpy [J/kg]
    inline scalar Ha(const scalar p, const scalar T) const;

    //- Sensible enthalpy [J/kg]
    inline scalar Hs(const scalar p, const scalar T) const;

    //- Chemical enthalpy [J/kg]
    inline scalar Hc() const;

    //- Absolute energy [J/kg]
    inline scalar Ea(const scalar p, const scalar T) const;

    //- Sensible energy [J/kg]
    inline scalar Es(const scalar p, const scalar T) const;

    //- Entropy [J/(kg K)]
    inline scalar S(const scalar p, const scalar T) const;


    // Derivative term used for Jacobian

    //- Derivative of Gibbs free energy w.r.t. temperature
    inline scalar dGdT(const scalar p, const scalar T) const;
    
    //- Temperature derivative of heat capacity at constant pressure
    inline scalar dCpdT(const scalar p, const scalar T) const;


    // I-O
    
    //- Write to Ostream
    void write(Ostream& os) const;
    
    
    // Member operators
    
    inline void operator+=(const CoolPropThermo&);
    
    
    // Friend operators
    
    friend CoolPropThermo operator+ <EquationOfState>
    (
        const CoolPropThermo&,
        const CoolPropThermo&
    );
    
    friend CoolPropThermo operator* <EquationOfState>
    (
        const scalar,
        const CoolPropThermo&
    );
    
    friend CoolPropThermo operator== <EquationOfState>
    (
        const CoolPropThermo&,
        const CoolPropThermo&
    );
    
    
    // Ostream Operator
    
    friend Ostream& operator<< <EquationOfState>
    (
        Ostream&,
        const CoolPropThermo&
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CoolPropThermoI.H"

#ifdef NoRepository
    #include "CoolPropThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
