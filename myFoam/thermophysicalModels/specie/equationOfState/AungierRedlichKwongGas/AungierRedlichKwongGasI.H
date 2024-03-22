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

#include "AungierRedlichKwongGas.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

#if (OPENFOAM >= 1912)
#define small SMALL
#endif

template<class Specie>
inline Foam::AungierRedlichKwongGas<Specie>::AungierRedlichKwongGas
(
    const Specie& sp,
    const scalar& Tc,
    const scalar& Vc,
    const scalar& Pc,
    const scalar& omega
)
:
    Specie(sp),
    Tc_(Tc),
    Vc_(Vc),
    Pc_(Pc),
    omega_(omega)
{
    n_ = 0.4986 + 1.1735*omega_ + 0.4754*sqr(omega_);
    a_ = 0.42747*sqr(this->R()*Tc_)/Pc_;
    b_ = 0.08664*this->R()*Tc_/Pc;
    c_ = this->R()*Tc_/(Pc_ + a_/(Vc_*(Vc_ + b_))) + b_ - Vc_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::AungierRedlichKwongGas<Specie>::AungierRedlichKwongGas
(
    const word& name,
    const AungierRedlichKwongGas& pg
)
:
    Specie(name, pg),
    Tc_(pg.Tc_),
    Vc_(pg.Vc_),
    Pc_(pg.Pc_),
    omega_(pg.omega_),
    a_(pg.a_),
    b_(pg.b_),
    c_(pg.c_),
    n_(pg.n_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::AungierRedlichKwongGas <Specie>>
Foam::AungierRedlichKwongGas<Specie>::clone() const
{
    return autoPtr<AungierRedlichKwongGas<Specie>>
    (
        new AungierRedlichKwongGas<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::AungierRedlichKwongGas<Specie>>
Foam::AungierRedlichKwongGas<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<AungierRedlichKwongGas<Specie>>
    (
        new AungierRedlichKwongGas<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::AungierRedlichKwongGas<Specie>::rho
(
    scalar p,
    scalar T
) const
{
    const scalar Z = this->Z(p, T);
    return p/(Z*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::AungierRedlichKwongGas<Specie>::dpdV(scalar V, scalar T, scalar alpha)
    const
{
    return (- this->R()*T/sqr(V - b_ + c_) + alpha/(V*sqr(V + b_)) + alpha/(sqr(V)*(V + b_)));
}


template<class Specie>
inline Foam::scalar Foam::AungierRedlichKwongGas<Specie>::dpdT(scalar V, scalar T, scalar alpha)
    const
{
    return this->R()/(V - b_ + c_) + alpha*n_/(T*V*(V + b_));
}


template<class Specie>
inline Foam::scalar Foam::AungierRedlichKwongGas<Specie>::H(scalar p, scalar T) const
{

    const scalar Z = this->Z(p, T);
    const scalar V = Z*this->R()*T/p;
    const scalar Tr = T/Tc_;

    scalar Hdep = p*V - this->R()*T - a_/b_*(n_ + 1)*pow(Tr, -n_)*log((V + b_)/V);
    return Hdep;
}


template<class Specie>
inline Foam::scalar Foam::AungierRedlichKwongGas<Specie>::Cp(scalar p, scalar T) const
{
    const scalar Z = this->Z(p, T);
    const scalar V = Z*this->R()*T/p;
    const scalar Tr = T/Tc_;
    const scalar alpha = a_*pow(Tr, -n_);

    const scalar dVdT = - this->dpdT(V, T, alpha)/this->dpdV(V, T, alpha);

    // Cp = dHdT_p
    const scalar cpdep = p*dVdT - this->R()
        + alpha/b_*n_*(n_ + 1)/T*log((V + b_)/V)
        + alpha*(n_ + 1)*dVdT/V/(V + b_);

    return cpdep;
}


template<class Specie>
inline Foam::scalar Foam::AungierRedlichKwongGas<Specie>::S
(
    scalar p,
    scalar T
) const
{
    const scalar Tr = T/Tc_;
    const scalar alpha = a_*pow(Tr, -n_);

    const scalar Z = this->Z(p, T);
    const scalar V = Z*this->R()*T/p;
    const scalar Vstd = this->R()*Tstd/Pstd;
    
    const scalar Sdep = -this->R()*log(T/Tstd) + this->R()*log((V - b_ + c_)/Vstd) - n_*alpha/b_/T*log((V + b_)/V);
    return Sdep;
}


template<class Specie>
inline Foam::scalar Foam::AungierRedlichKwongGas<Specie>::psi
(
    scalar p,
    scalar T
) const
{
    const scalar Z = this->Z(p, T);
    return 1.0/(Z*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::AungierRedlichKwongGas<Specie>::Z
(
    scalar p,
    scalar T
) const
{
    const scalar Tr = T/Tc_;
    const scalar RT = this->R()*T;
    const scalar alpha = a_*pow(Tr, -n_);
    
    const scalar B = c_*p/RT - 1;
    const scalar C = p*(-RT*b_ + alpha - sqr(b_)*p + b_*c_*p)/sqr(RT);
    const scalar D = alpha*sqr(p)*(c_ - b_)/pow3(RT);

    // Solve Z^3 + B*Z^2 + C*Z + D = 0, return maximal real root 
    scalar root = 0;
    
    const scalar q = (3*C - B*B)/9;
    const scalar r = (-27*D + B*(9*C - 2*B*B))/54.0;
    const scalar disc = pow3(q) + sqr(r);
    const scalar term1 = B/3.0;
    
    if (disc > 0)
    {   // one root real, two are complex
        scalar s = cbrt(r + sqrt(disc));
        scalar t = cbrt(r - sqrt(disc));
        root = -term1 + s + t;
    } 
    else if (disc == 0)
    {   // All roots real, at least two are equal.
        scalar r13 = std::cbrt(r); //((r < 0) ? -pow(-r,1.0/3.0) : pow(r,1.0/3.0));
        root = max(2*r13 - term1, -r13 - term1);
    }
    else
    {
        // Only option left is that all roots are real and unequal (to get here, q < 0)
        const scalar dum1 = acos(r/sqrt(-pow3(q)));
        const scalar r13 = 2.0*sqrt(-q);
        const scalar x1 = -term1 + r13*cos(dum1/3.0);
        const scalar x2 = -term1 + r13*cos((dum1 + 2.0*constant::mathematical::pi)/3.0);
        const scalar x3 = -term1 + r13*cos((dum1 + 4.0*constant::mathematical::pi)/3.0);

        root = max(x1, max(x2, x3));
    }

    return root;
}


template<class Specie>
inline Foam::scalar Foam::AungierRedlichKwongGas<Specie>::CpMCv
(
    scalar p,
    scalar T
) const
{
    const scalar Z = this->Z(p, T);
    const scalar V = Z*this->R()*T/p;
    const scalar Tr = T/Tc_;
    const scalar alpha = a_*pow(Tr, -n_);

    const scalar kappa= -1/(V*this->dpdV(V,T,alpha));
    const scalar beta = this->dpdT(V, T, alpha)*kappa;

    return V*T*sqr(beta)/kappa;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::AungierRedlichKwongGas<Specie>::operator+=
(
    const AungierRedlichKwongGas<Specie>& pg
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
inline void Foam::AungierRedlichKwongGas<Specie>::operator*=(const scalar s)
{
     Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


template<class Specie>
Foam::AungierRedlichKwongGas<Specie> Foam::operator+
(
    const AungierRedlichKwongGas<Specie>& pg1,
    const AungierRedlichKwongGas<Specie>& pg2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pg1)
      + static_cast<const Specie&>(pg2)
    );

    if (mag(sp.Y()) < small)
    {
        return AungierRedlichKwongGas<Specie>
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
        
        return AungierRedlichKwongGas<Specie>
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
Foam::AungierRedlichKwongGas<Specie> Foam::operator*
(
    const scalar s,
    const AungierRedlichKwongGas<Specie>& pg
)
{
    return AungierRedlichKwongGas<Specie>
    (
        s*static_cast<const Specie&>(pg),
        pg.Tc_,
        pg.Vc_,
        pg.Pc_,
        pg.omega_
    );
}


template<class Specie>
Foam::AungierRedlichKwongGas<Specie> Foam::operator==
(
    const AungierRedlichKwongGas<Specie>& pg1,
    const AungierRedlichKwongGas<Specie>& pg2
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

    return AungierRedlichKwongGas<Specie>
    (
        sp,
        Tc,
        Vc,
        Pc,
        Y2*pg2.omega_ - Y1*pg1.omega_
    );
}


// ************************************************************************* //
