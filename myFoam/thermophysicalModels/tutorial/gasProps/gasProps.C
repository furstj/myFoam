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

Application
    gasProps

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "gasProperties.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    autoPtr<fluidThermo> pthermo(fluidThermo::New(mesh));

    autoPtr<gasProperties> gp(gasProperties::New(pthermo()));

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    scalar p = 101325;
    scalar T = 288.15;
    
    Info << "W          = " << gp->W() << endl;
    Info << "rho(p,T)   = " << gp->rho(p,T) << endl;
    Info << "psi(p,T)   = " << gp->psi(p,T) << endl;
    Info << "CpMCv(p,T) = " << gp->CpMCv(p,T) << endl;
    Info << "Cp(p,T)    = " << gp->Cp(p,T) << endl;
    Info << "Ha(p,T)    = " << gp->Ha(p,T) << endl;
    Info << "Hs(p,T)    = " << gp->Hs(p,T) << endl;
    Info << "Hc()       = " << gp->Hc() << endl;
    Info << "S(p,T)     = " << gp->S(p,T) << endl;
    Info << "mu(p,T)    = " << gp->mu(p,T) << endl;
    Info << "kappa(p,T) = " << gp->kappa(p,T) << endl;
    Info << "alphah(p,T)= " << gp->alphah(p,T) << endl;
    
    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
