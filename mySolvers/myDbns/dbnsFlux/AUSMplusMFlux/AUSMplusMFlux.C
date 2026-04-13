/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "AUSMplusMFlux.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(AUSMplusMFlux, 0);
    addToRunTimeSelectionTable(dbnsFlux, AUSMplusMFlux, dictionary);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::AUSMplusMFlux::AUSMplusMFlux(const fvMesh&, const dictionary& dict)
{
    dictionary    mySubDict( dict.subOrEmptyDict("AUSMplusMFluxCoeffs") );
    sqrMachInf_ = mySubDict.lookupOrDefault("sqrMachInf", 0.01);
    alpha0_     = mySubDict.lookupOrDefault("alpha0", 3.0/16.0);
    requiresPressureSensing_ = true;
    
    if (mySubDict.lookupOrDefault("printCoeffs", false))
        Info << mySubDict << nl;
};


void Foam::AUSMplusMFlux::evaluateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const scalar& pLeft,
    const scalar& pRight,
    const vector& ULeft,
    const vector& URight,
    const scalar& TLeft,
    const scalar& TRight,
    const scalar& RLeft,
    const scalar& RRight,
    const scalar& CvLeft,
    const scalar& CvRight,
    const vector& Sf,
    const scalar& magSf,
    const scalar& meshPhi,
    const scalar& h
) const
{
    // Step 1: decode left and right:
    // normal vector
    const vector normalVector = Sf/magSf;

    // Ratio of specific heat capacities
    const scalar kappaLeft = (RLeft + CvLeft)/CvLeft;
    const scalar kappaRight = (RRight + CvRight)/CvRight;
    const scalar kappa = 0.5 * (kappaLeft + kappaRight);

    // Density
    const scalar rhoLeft = pLeft/(RLeft*TLeft);
    const scalar rhoRight = pRight/(RRight*TRight);
    const scalar rho_mean = 0.5*(rhoLeft + rhoRight);
    // DensityTotalEnthalpy
    const scalar rhoHLeft  = rhoLeft*(CvLeft*TLeft + 0.5*magSqr(ULeft)) + pLeft;
    const scalar rhoHRight = rhoRight*(CvRight*TRight + 0.5*magSqr(URight)) + pRight;

    const scalar u_normal_L = (ULeft & normalVector);
    const scalar u_normal_R = (URight & normalVector);
    
    const scalar u_tangential_L = mag(ULeft - (ULeft & normalVector) * normalVector);
    const scalar u_tangential_R = mag(URight - (URight & normalVector) * normalVector);

    
    const scalar h_L = rhoHLeft/rhoLeft;
    const scalar h_R = rhoHRight/rhoRight;
    
    // calculate normal component of enthalpy
    const scalar h_normal = 0.5 * (h_L - 0.5 * sqr(u_tangential_L) + h_R - 0.5 * sqr(u_tangential_R));


    const scalar sqr_c_star = 2.0 * (kappa - 1.0) / (kappa + 1.0) * h_normal;
    
    const scalar u_normal_mean = 0.5 * (u_normal_L + u_normal_R);

    const scalar c_face = pos0(u_normal_mean) * (sqr_c_star / max(mag(u_normal_L), sqrt(sqr_c_star))) + neg(u_normal_mean) * (sqr_c_star / max(mag(u_normal_R), sqrt(sqr_c_star)));

    const scalar Mach_L = u_normal_L / c_face;
    const scalar Mach_R = u_normal_R / c_face;
    const scalar magnitude_Mach_L = mag(Mach_L);
    const scalar magnitude_Mach_R = mag(Mach_R);
    
    const scalar mach = min(1.0, max(magnitude_Mach_L, magnitude_Mach_R));
    
    const scalar PI_ = constant::mathematical::pi;
    const scalar f = 0.5 * (1 - cos(PI_ * mach));

    const scalar f_0 = min(1.0, max(f, sqrMachInf_));

    const scalar alpha = alpha0_;

    const scalar g = 0.5 * (1.0 + cos(PI_ * h));

    const scalar p_plus_L = pos0(magnitude_Mach_L - 1.0) * 0.5 * (1.0 + sign(Mach_L)) + neg(magnitude_Mach_L - 1.0) * (0.25 * sqr(Mach_L + 1.0) * (2.0 - Mach_L) + alpha * Mach_L * sqr(sqr(Mach_L) - 1.0));
    
    const scalar p_minus_R = pos0(magnitude_Mach_R - 1.0) * 0.5 * (1.0 - sign(Mach_R)) + neg(magnitude_Mach_R - 1.0) * (0.25 * sqr(Mach_R - 1.0) * (2.0 + Mach_R) - alpha * Mach_R * sqr(sqr(Mach_R) - 1.0));
    
    const scalar p_mean = 0.5 * (pLeft + pRight);
    const scalar dp = (p_plus_L - p_minus_R) * 0.5 * (pLeft - pRight) + f_0 * ((p_plus_L + p_minus_R - 1) * p_mean);

    const scalar ps = p_mean + dp;

    const vector pu = -g * kappa * p_mean / c_face * p_plus_L * p_minus_R * (URight - ULeft);

    const vector p12 = ps * normalVector + pu;

    const scalar Mach_plus_L = pos0(magnitude_Mach_L - 1.0) * 0.5 * (Mach_L + magnitude_Mach_L) + neg(magnitude_Mach_L - 1.0) * (0.25 * sqr(Mach_L + 1.0) + 0.125 * sqr(sqr(Mach_L) - 1.0));
    const scalar Mach_minus_R = pos0(magnitude_Mach_R - 1.0) * 0.5 * (Mach_R - magnitude_Mach_R) + neg(magnitude_Mach_R - 1.0) * (-0.25 * sqr(Mach_R - 1.0) - 0.125 * sqr(sqr(Mach_R) - 1.0));

    const scalar Mp = -0.5 * (1.0 - f) * (pRight - pLeft) / rho_mean / sqr(c_face) * (1.0 - g);

    const scalar Mach_1_2 = Mach_plus_L + Mach_minus_R + Mp;

    const scalar mdot = c_face * Mach_1_2 * (rhoLeft * pos0(Mach_1_2) + rhoRight * neg(Mach_1_2));

    rhoFlux  = magSf * mdot;
    rhoUFlux = magSf * (0.5 * (mdot + mag(mdot)) * ULeft + 0.5 * (mdot - mag(mdot)) * URight + p12);
    rhoEFlux = magSf * (0.5 * (mdot + mag(mdot)) * h_L + 0.5 * (mdot - mag(mdot)) * h_R);
}


// ************************************************************************* //
