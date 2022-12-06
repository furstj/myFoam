/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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

namespace Foam
{


template<typename TransportModel>
genericGasProperties<TransportModel>::genericGasProperties(const dictionary& dict):
    //gasProperties(dict),
    transport_(dict.subDict("mixture"))
{
}


template<typename TransportModel>
scalar genericGasProperties<TransportModel>::W() const
{
    return transport_.W();
}


template<typename TransportModel>
scalar genericGasProperties<TransportModel>::R() const
{
    return transport_.R();
}


template<typename TransportModel>
scalar genericGasProperties<TransportModel>::rho(scalar p, scalar T) const
{
    return transport_.rho(p,T);
}


template<typename TransportModel>
scalar genericGasProperties<TransportModel>::Z(scalar p, scalar T) const
{
    return transport_.Z(p,T);
}


template<typename TransportModel>
scalar genericGasProperties<TransportModel>::psi(scalar p, scalar T) const
{
    return transport_.psi(p,T);
}


template<typename TransportModel>
scalar genericGasProperties<TransportModel>::CpMCv(scalar p, scalar T) const
{
    return transport_.CpMCv(p,T);
}


template<typename TransportModel>
scalar genericGasProperties<TransportModel>::Cp(scalar p, scalar T) const
{
    return transport_.Cp(p,T);
}


template<typename TransportModel>
scalar genericGasProperties<TransportModel>::Cv(scalar p, scalar T) const
{
    return transport_.Cv(p,T);
}


template<typename TransportModel>
scalar genericGasProperties<TransportModel>::Ha(scalar p, scalar T) const
{
    return transport_.Ha(p,T);
}


template<typename TransportModel>
scalar genericGasProperties<TransportModel>::Hs(scalar p, scalar T) const
{
    return transport_.Hs(p,T);
}


template<typename TransportModel>
scalar genericGasProperties<TransportModel>::Hc() const
{
    return transport_.Hc();
}


template<typename TransportModel>
scalar genericGasProperties<TransportModel>::Ea(scalar p, scalar T) const
{
    return transport_.Ea(p,T);
}


template<typename TransportModel>
scalar genericGasProperties<TransportModel>::Es(scalar p, scalar T) const
{
    return transport_.Es(p,T);
}


template<typename TransportModel>
scalar genericGasProperties<TransportModel>::S(scalar p, scalar T) const
{
    return transport_.S(p,T);
}

template<typename TransportModel>
scalar genericGasProperties<TransportModel>::beta_p(scalar p, scalar T) const
{
    scalar v  = 1/rho(p,T);
    scalar T1 = (1 + 1.e-6)*T;
    scalar v1 = 1/rho(p,T1);
    return 1/v*(v1 - v)/(T1 - T);
}

template<typename TransportModel>
scalar genericGasProperties<TransportModel>::beta_T(scalar p, scalar T) const
{
    scalar v  = 1/rho(p,T);
    scalar p1 = (1 + 1e-6)*p;
    scalar v1 = 1/rho(p1,T);    
    return -1/v*(v1 - v)/(p1 - p);
}

template<typename TransportModel>
scalar genericGasProperties<TransportModel>::c(scalar p, scalar T) const
{
    scalar cp = Cp(p,T);
    scalar cv = cp - CpMCv(p,T);
    scalar gamma_pv = cp/cv/beta_T(p,T)/p;
    scalar z = Z(p,T);
    return sqrt(gamma_pv*z*R()*T);
}


template<typename TransportModel>
scalar genericGasProperties<TransportModel>::THs(const scalar Hs, const scalar p, const scalar T0)
    const
{
    return transport_.THs(Hs, p, T0);
}
    
template<typename TransportModel>
scalar genericGasProperties<TransportModel>::THa(const scalar Ha, const scalar p, const scalar T0)
    const
{
    return transport_.THa(Ha, p, T0);
}
    

template<typename TransportModel>
scalar genericGasProperties<TransportModel>::TEs(const scalar Es, const scalar p, const scalar T0)
    const
{
    return transport_.TEs(Es, p, T0);
}
    
template<typename TransportModel>
scalar genericGasProperties<TransportModel>::TEa(const scalar Ea, const scalar p, const scalar T0)
    const
{
    return transport_.TEa(Ea, p, T0);
}

template<typename TransportModel>
scalar genericGasProperties<TransportModel>::TpS(const scalar p, const scalar S, const scalar T0)
    const
{
    scalar T = T0;
    const scalar tol = 1.e-8;
    const label maxIter = 100;

    const scalar TTol = T*tol;
    
    label iter = 0;
    scalar dT;
    
    do
    {
        scalar dS = S - this->S(p, T);
        dT = T/this->Cp(p,T)*dS;
        T += dT;

        if (iter++ > maxIter)
        {
            FatalErrorInFunction
                << "Maximum number of iterations exceeded: " << maxIter
                << " when starting from T0:" << T0
                << " T  : " << T
                << " p  : " << p
                << " tol: " << tol
                << abort(FatalError);
        }
    } while ( mag(dT) > TTol);
    return T;
}


template<typename TransportModel>
scalar genericGasProperties<TransportModel>::pEs(const scalar Es, const scalar rho, const scalar p0)
    const
{
    scalar p = p0;
    scalar T = p0/(rho*R());   // Estimate from perfect gas!
    const scalar tol = 1.e-8;
    const label maxIter = 100;

    const scalar TTol = T*tol;
    const scalar pTol = p*tol;
    
    label iter = 0;
    scalar dT, dp;
    
    do
    {
        dT = (Es - this->Es(p, T))/ this->Cv(p,T);
        dp = this->Z(p,T)*R()*rho*dT;
        T += dT;
        p += dp;

        if (iter++ > maxIter)
        {
            FatalErrorInFunction
                << "Maximum number of iterations exceeded: " << maxIter
                << " when starting from p0:" << p0
                << " T  : " << T
                << " p  : " << p
                << " tol: " << tol
                << abort(FatalError);
        }
    } while ( (mag(dT) > TTol) || (mag(dp) > pTol) );
    return p;
}


template<typename TransportModel>
scalar genericGasProperties<TransportModel>::pHS(const scalar Hs, const scalar S, const scalar p0)
    const
{
    const scalar tol = 1.e-8;
    const label maxIter = 100;

    const scalar pTol = p0*tol;
    scalar T = THs(Hs, p0, 288.15);
    
    label iter = 0;
    scalar dp;
    scalar p = p0;
    do
    {
        T = TpS(p, S, T);
        scalar dh = Hs - this->Hs(p, T);
        scalar rho = this->rho(p, T);
        dp = rho*dh;
        p += dp;

        if (iter++ > maxIter)
        {
            FatalErrorInFunction
                << "Maximum number of iterations exceeded: " << maxIter
                << " when starting from p0:" << p0
                << " T  : " << T
                << " p  : " << p
                << " tol: " << tol
                << abort(FatalError);
        }
    } while ( mag(dp) > pTol);
    return p;
}


template<typename TransportModel>
scalar genericGasProperties<TransportModel>::mu(scalar p, scalar T) const
{
    return transport_.mu(p,T);
}


template<typename TransportModel>
scalar genericGasProperties<TransportModel>::kappa(scalar p, scalar T) const
{
    return transport_.kappa(p,T);
}


template<typename TransportModel>
scalar genericGasProperties<TransportModel>::alphah(scalar p, scalar T) const
{
    return transport_.alphah(p,T);
}

template<typename TransportModel>
void genericGasProperties<TransportModel>::writeData(Ostream& os) const
{}

}



// ************************************************************************* //
