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

\*---------------------------------------------------------------------------*/

#include "pVirialGas.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

#if (OPENFOAM >= 1912)
#define small SMALL
#endif

template<class Specie>
inline Foam::pVirialGas<Specie>::pVirialGas
(
    const Specie& sp,
    const scalar& z0,
    const scalar& z1
)
:
    Specie(sp),
    z0_(z0),
    z1_(z1)
{
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::pVirialGas<Specie>::pVirialGas
(
    const word& name,
    const pVirialGas& pg
)
:
    Specie(name, pg),
    z0_(pg.z0_),
    z1_(pg.z1_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::pVirialGas <Specie>>
Foam::pVirialGas<Specie>::clone() const
{
    return autoPtr<pVirialGas<Specie>>
    (
        new pVirialGas<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::pVirialGas<Specie>>
Foam::pVirialGas<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<pVirialGas<Specie>>
    (
        new pVirialGas<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::pVirialGas<Specie>::p
(
    scalar rho,
    scalar T
) const
{
    const scalar rhoRT = rho*this->R()*T;
    return rhoRT*z0_ / (1 - z1_*rhoRT);
}

template<class Specie>
inline Foam::scalar Foam::pVirialGas<Specie>::rho
(
    scalar p,
    scalar T
) const
{
    const scalar Z = this->Z(p, T);
    return p/(Z*this->R()*T);
}



template<class Specie>
inline Foam::scalar Foam::pVirialGas<Specie>::H(scalar p, scalar T) const
{
    scalar Hdep = 0.0;
    return Hdep;
}


template<class Specie>
inline Foam::scalar Foam::pVirialGas<Specie>::Cp(scalar p, scalar T) const
{
    scalar cpdep = 0.0;
    return cpdep;
}


template<class Specie>
inline Foam::scalar Foam::pVirialGas<Specie>::S
(
    scalar p,
    scalar T
) const
{
    const scalar Sdep = -this->R()*(z0_*log(p/Pstd) + z1_*p);
    return Sdep;
}


template<class Specie>
inline Foam::scalar Foam::pVirialGas<Specie>::psi
(
    scalar p,
    scalar T
) const
{
    const scalar Z = this->Z(p, T);
    return 1.0/(Z*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::pVirialGas<Specie>::Z
(
    scalar p,
    scalar T
) const
{
    const scalar z = z0_ + z1_*p;
    return z;
}


template<class Specie>
inline Foam::scalar Foam::pVirialGas<Specie>::CpMCv
(
    scalar p,
    scalar T
) const
{
    return this->R()*sqr(z0_ + z1_*p)/z0_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::pVirialGas<Specie>::operator+=
(
    const pVirialGas<Specie>& pg
)
{
    scalar Y1 = this->Y();
    Specie::operator+=(pg);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        const scalar Y2 = pg.Y()/this->Y();

        z0_ = Y1*z0_ + Y2*pg.z0_;
        z1_ = Y1*z1_ + Y2*pg.z1_;
    }
}


template<class Specie>
inline void Foam::pVirialGas<Specie>::operator*=(const scalar s)
{
     Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


template<class Specie>
Foam::pVirialGas<Specie> Foam::operator+
(
    const pVirialGas<Specie>& pg1,
    const pVirialGas<Specie>& pg2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pg1)
      + static_cast<const Specie&>(pg2)
    );

    if (mag(sp.Y()) < small)
    {
        return pVirialGas<Specie>
        (
            sp,
            pg1.z0_,
            pg1.z1_
        );
    }
    else
    {
        const scalar Y1 = pg1.Y()/sp.Y();
        const scalar Y2 = pg2.Y()/sp.Y();
        
        return pVirialGas<Specie>
        (
            sp,
            Y1*pg1.z0_ + Y2*pg2.z0_,
            Y1*pg1.z1_ + Y2*pg2.z1_
        );
    }
}


template<class Specie>
Foam::pVirialGas<Specie> Foam::operator*
(
    const scalar s,
    const pVirialGas<Specie>& pg
)
{
    return pVirialGas<Specie>
    (
        s*static_cast<const Specie&>(pg),
        pg.z0_,
        pg.z1_
    );
}


template<class Specie>
Foam::pVirialGas<Specie> Foam::operator==
(
    const pVirialGas<Specie>& pg1,
    const pVirialGas<Specie>& pg2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pg1)
     == static_cast<const Specie&>(pg2)
    );

    const scalar Y1 = pg1.Y()/sp.Y();
    const scalar Y2 = pg2.Y()/sp.Y();

    const scalar z0 = Y2*pg2.z0_ - Y1*pg1.z0_;
    const scalar z1 = Y2*pg2.z1_ - Y1*pg1.z1_;

    return pVirialGas<Specie>
    (
        sp,
        z0,
        z1
    );
}


// ************************************************************************* //
