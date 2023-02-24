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

#include "CoolPropGas.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


template<class Specie>
inline Foam::CoolPropGas<Specie>::CoolPropGas
(
    const Specie& sp,
    const std::string& backend,
    const std::string& fluid_names
)
:
    Specie(sp),
    state_(CoolProp::AbstractState::factory(backend, fluid_names))
{
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class Specie>
inline Foam::CoolPropGas<Specie>::CoolPropGas
(
    const word& name,
    const CoolPropGas& pg
)
:
    Specie(name, pg),
    state_(CoolProp::AbstractState::factory
        (pg.state_->backend_name(), pg.state_->name()))
{}


template<class Specie>
inline Foam::autoPtr<Foam::CoolPropGas <Specie>>
Foam::CoolPropGas<Specie>::clone() const
{
    return autoPtr<CoolPropGas<Specie>>
    (
        new CoolPropGas<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::CoolPropGas<Specie>>
Foam::CoolPropGas<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<CoolPropGas<Specie>>
    (
        new CoolPropGas<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::CoolPropGas<Specie>::rho
(
    scalar p,
    scalar T
) const
{
    this->state_->update(CoolProp::PT_INPUTS, p, T);
    return this->state_->rhomass();
}


template<class Specie>
inline Foam::scalar Foam::CoolPropGas<Specie>::H(scalar p, scalar T) const
{
    this->state_->update(CoolProp::PT_INPUTS, p, T);
    return this->state_->hmass();
}

template<class Specie>
inline Foam::scalar Foam::CoolPropGas<Specie>::E(scalar p, scalar T) const
{
    this->state_->update(CoolProp::PT_INPUTS, p, T);
    return this->state_->umass();
}

template<class Specie>
inline Foam::scalar Foam::CoolPropGas<Specie>::Cp(scalar p, scalar T) const
{
    this->state_->update(CoolProp::PT_INPUTS, p, T);
    return this->state_->cpmass();
}

template<class Specie>
inline Foam::scalar Foam::CoolPropGas<Specie>::Cv(scalar p, scalar T) const
{
    this->state_->update(CoolProp::PT_INPUTS, p, T);
    return this->state_->cvmass();
}

template<class Specie>
inline Foam::scalar Foam::CoolPropGas<Specie>::S(scalar p, scalar T) const
{
    this->state_->update(CoolProp::PT_INPUTS, p, T);
    return this->state_->smass();
 }
 

template<class Specie>
inline Foam::scalar Foam::CoolPropGas<Specie>::psi(scalar p, scalar T) const
{
    scalar rho = this->rho(p,T);
    return rho/p;
}


template<class Specie>
inline Foam::scalar Foam::CoolPropGas<Specie>::Z(scalar p, scalar T) const
{
    scalar rho = this->rho(p,T);
    return p/(rho*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::CoolPropGas<Specie>::CpMCv(scalar p, scalar T) const
{
    this->state_->update(CoolProp::PT_INPUTS, p, T);
    return this->state_->cpmass() - this->state_->cvmass();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

/*
template<class Specie>
inline void Foam::CoolPropGas<Specie>::operator+=
(
    const CoolPropGas<Specie>& pg
)
{
    scalar Y1 = this->Y();
    Specie::operator+=(pg);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        const scalar Y2 = pg.Y()/this->Y();

        Tc_ = Y1*Tc_ + Y2*pg.Tc_;
        Vc_ = Y1*Vc_ + Y2*pg.Vc_;
        Pc_ = Y1*Pc_ + Y2*pg.Pc_;
        omega_ = Y1*omega_ + Y2*pg.omega_;
        n_ = 0.4986 + 1.1735*omega_ + 0.4754*sqr(omega_);
        a_ = 0.42747*sqr(this->R()*Tc_)/Pc_;
        b_ = 0.08664*this->R()*Tc_/Pc_;
        c_ = this->R()*Tc_/(Pc_ + a_/(Vc_*(Vc_ + b_))) + b_ - Vc_;
    }
}


template<class Specie>
inline void Foam::CoolPropGas<Specie>::operator*=(const scalar s)
{
     Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


template<class Specie>
Foam::CoolPropGas<Specie> Foam::operator+
(
    const CoolPropGas<Specie>& pg1,
    const CoolPropGas<Specie>& pg2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pg1)
      + static_cast<const Specie&>(pg2)
    );

    if (mag(sp.Y()) < small)
    {
        return CoolPropGas<Specie>
        (
            sp,
            pg1.Tc_,
            pg1.Vc_,
            pg1.Pc_,
            pg1.omega_
        );
    }
    else
    {
        const scalar Y1 = pg1.Y()/sp.Y();
        const scalar Y2 = pg2.Y()/sp.Y();

        const scalar Z1 = pg1.Pc_*pg1.Vc_/(pg1.R()*pg1.Tc_);
        const scalar Z2 = pg2.Pc_*pg2.Vc_/(pg2.R()*pg2.Tc_);
    
        const scalar Tc = Y1*pg1.Tc_ + Y2*pg2.Tc_;
        const scalar Vc = (Y1*pg1.Vc_*pg1.W() + Y2*pg2.Vc_*pg2.W())/sp.W();
        //const scalar Pc = Y1*pg1.Pc_ + Y2*pg2.Pc_;
        const scalar Zc = Y1*Z1 + Y2*Z2;
        
        return CoolPropGas<Specie>
        (
            sp,
            Tc,
            Vc,
            Zc*sp.R()*Tc/Vc,
            Y1*pg1.omega_ + Y2*pg2.omega_
        );
    }
}


template<class Specie>
Foam::CoolPropGas<Specie> Foam::operator*
(
    const scalar s,
    const CoolPropGas<Specie>& pg
)
{
    return CoolPropGas<Specie>
    (
        s*static_cast<const Specie&>(pg),
        pg.Tc_,
        pg.Vc_,
        pg.Pc_,
        pg.omega_
    );
}


template<class Specie>
Foam::CoolPropGas<Specie> Foam::operator==
(
    const CoolPropGas<Specie>& pg1,
    const CoolPropGas<Specie>& pg2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pg1)
     == static_cast<const Specie&>(pg2)
    );

    const scalar Y1 = pg1.Y()/sp.Y();
    const scalar Y2 = pg2.Y()/sp.Y();

    const scalar Tc = Y2*pg2.Tc_ - Y1*pg1.Tc_;
    const scalar Vc = Y2*pg2.Vc_ - Y1*pg1.Vc_;
    const scalar Pc = Y2*pg2.Pc_ - Y1*pg1.Pc_;

    return CoolPropGas<Specie>
    (
        sp,
        Tc,
        Vc,
        Pc,
        Y2*pg2.omega_ - Y1*pg1.omega_
    );
}
*/

// ************************************************************************* //
