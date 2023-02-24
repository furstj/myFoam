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

#include "specie.H"
#include "CoolProp.h"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::CoolPropTransport<Thermo>::CoolPropTransport
(
    const Thermo& t
)
:
    Thermo(t)
{}


template<class Thermo>
inline Foam::CoolPropTransport<Thermo>::CoolPropTransport
(
    const word& name,
    const CoolPropTransport& st
)
:
    Thermo(name, st)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::CoolPropTransport<Thermo>>
Foam::CoolPropTransport<Thermo>::clone() const
{
    return autoPtr<CoolPropTransport<Thermo>>
    (
        new CoolPropTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::CoolPropTransport<Thermo>>
Foam::CoolPropTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<CoolPropTransport<Thermo>>
    (
        new CoolPropTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::CoolPropTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    auto state = Thermo::CoolPropGas::state();
    state->update(CoolProp::PT_INPUTS, p, T);
    return state->viscosity();
}


template<class Thermo>
inline Foam::scalar Foam::CoolPropTransport<Thermo>::kappa
(
    const scalar p, const scalar T
) const
{
    auto state = Thermo::CoolPropGas::state();
    state->update(CoolProp::PT_INPUTS, p, T);
    return state->conductivity();
}


template<class Thermo>
inline Foam::scalar Foam::CoolPropTransport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{
    auto state = Thermo::CoolPropGas::state();
    state->update(CoolProp::PT_INPUTS, p, T);
    return state->conductivity()/state->cpmass();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::CoolPropTransport<Thermo>::operator=
(
    const CoolPropTransport<Thermo>& st
)
{
    Thermo::operator=(st);
}


template<class Thermo>
inline void Foam::CoolPropTransport<Thermo>::operator+=
(
    const CoolPropTransport<Thermo>& st
)
{
    Thermo::operator+=(st);
}


template<class Thermo>
inline void Foam::CoolPropTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::CoolPropTransport<Thermo> Foam::operator+
(
    const CoolPropTransport<Thermo>& st1,
    const CoolPropTransport<Thermo>& st2
)
{
    return CoolPropTransport<Thermo>
    (
        static_cast<const Thermo&>(st1) + static_cast<const Thermo&>(st2)
    );
}


template<class Thermo>
inline Foam::CoolPropTransport<Thermo> Foam::operator*
(
    const scalar s,
    const CoolPropTransport<Thermo>& st
)
{
    return CoolPropTransport<Thermo>
    (
        s*static_cast<const Thermo&>(st)
    );
}


// ************************************************************************* //
