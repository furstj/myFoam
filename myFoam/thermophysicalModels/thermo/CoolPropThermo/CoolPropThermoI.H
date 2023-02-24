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

\*---------------------------------------------------------------------------*/

#include "CoolPropThermo.H"
#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::CoolPropThermo<EquationOfState>::CoolPropThermo
(
    const EquationOfState& st
)
:
    EquationOfState(st)
{}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::CoolPropThermo<EquationOfState>::CoolPropThermo
(
    const word& name,
    const CoolPropThermo& jt
)
:
    EquationOfState(name, jt)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::CoolPropThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::CoolPropThermo<EquationOfState>::Cp
(
    const scalar p,
    const scalar T
) const
{
    return EquationOfState::Cp(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::CoolPropThermo<EquationOfState>::Cv
(
    const scalar p,
    const scalar T
) const
{
    return EquationOfState::Cv(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::CoolPropThermo<EquationOfState>::Ha
(
    const scalar p,
    const scalar T
) const
{
    return EquationOfState::H(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::CoolPropThermo<EquationOfState>::Hs
(
    const scalar p,
    const scalar T
) const
{
    return Ha(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::CoolPropThermo<EquationOfState>::Hc() const
{
    return 0;
}

template<class EquationOfState>
inline Foam::scalar Foam::CoolPropThermo<EquationOfState>::Ea
(
    const scalar p,
    const scalar T
) const
{
    return EquationOfState::E(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::CoolPropThermo<EquationOfState>::Es
(
    const scalar p,
    const scalar T
) const
{
    return Ea(p, T);
}

template<class EquationOfState>
inline Foam::scalar Foam::CoolPropThermo<EquationOfState>::S
(
    const scalar p,
    const scalar T
) const
{
    return EquationOfState::S(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::CoolPropThermo<EquationOfState>::dGdT
(
    const scalar p,
    const scalar T
) const
{
    return EquationOfState::dGdT;
}


template<class EquationOfState>
inline Foam::scalar Foam::CoolPropThermo<EquationOfState>::dCpdT
(
    const scalar p,
    const scalar T
) const
{
    return EquationOfState::dcpdT(p,T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::CoolPropThermo<EquationOfState>::operator+=
(
    const CoolPropThermo<EquationOfState>& jt
)
{
    EquationOfState::operator+=(jt);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::CoolPropThermo<EquationOfState> Foam::operator+
(
    const CoolPropThermo<EquationOfState>& jt1,
    const CoolPropThermo<EquationOfState>& jt2
)
{
    EquationOfState eofs = jt1;
    eofs += jt2;
    
    return CoolPropThermo<EquationOfState>(eofs);
}


template<class EquationOfState>
inline Foam::CoolPropThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const CoolPropThermo<EquationOfState>& jt
)
{
    return CoolPropThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(jt)
    );
}


template<class EquationOfState>
inline Foam::CoolPropThermo<EquationOfState> Foam::operator==
(
    const CoolPropThermo<EquationOfState>& jt1,
    const CoolPropThermo<EquationOfState>& jt2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(jt1)
     == static_cast<const EquationOfState&>(jt2)
    );

    return CoolPropThermo<EquationOfState>
    (
        eofs
    );
}


// ************************************************************************* //
