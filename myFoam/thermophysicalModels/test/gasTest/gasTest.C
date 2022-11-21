/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 AUTHOR,AFFILIATION
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
    gasTest

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "gasProperties.H"
#include "fluidThermo.H"
#include <fstream>

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

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Range for p-T table
    scalar p1 = 1e3;
    label  np = 50;
    scalar kp = 1.15;

    scalar T1 = 200;
    label  nT = 40;
    scalar dT = 10;

    auto out = std::ofstream("thermo.vtk");
    out << "# vtk DataFile Version 2.0" << nl;
    out << "p-T diagrams" << nl;
    out << "ASCII" << nl;
    out << "DATASET RECTILINEAR_GRID" << nl;
    out << "DIMENSIONS " << nT << " " << np << " 1" << nl;
    out << "X_COORDINATES " << nT << " double" << nl;
    scalar T = T1;
    for (label iT = 0; iT < nT; iT++)
    {
        out << T << " ";
        T += dT;
    }
    out << nl;
    out << "Y_COORDINATES " << np << " double" << nl;
    scalar p = p1;
    for (label ip = 0; ip < np; ip++)
    {
        out << p/1000 << " ";
        p *= kp;
    }
    out << nl;

    out << "Z_COORDINATES " << 1 << " double" << nl;
    out << 0 << nl;

    out << "POINT_DATA " << np*nT << nl;

    out << "SCALARS density double 1" << nl;
    out << "LOOKUP_TABLE default" << nl;
    p = p1;
    for (label ip = 0; ip < np; ip++)
    {
        T = T1;
        for (label iT = 0; iT < nT; iT++)
        {
            out << gp->rho(p,T) << nl;
            T += dT;
        }
        p *= kp;
    }

    out << "SCALARS enthalpy double 1" << nl;
    out << "LOOKUP_TABLE default" << nl;
    p = p1;
    for (label ip = 0; ip < np; ip++)
    {
        T = T1;
        for (label iT = 0; iT < nT; iT++)
        {
            out << gp->Hs(p,T) << nl;
            T += dT;
        }
        p *= kp;
    }

    out << "SCALARS entropy double 1" << nl;
    out << "LOOKUP_TABLE default" << nl;
    p = p1;
    for (label ip = 0; ip < np; ip++)
    {
        T = T1;
        for (label iT = 0; iT < nT; iT++)
        {
            out << gp->S(p,T) << nl;
            T += dT;
        }
        p *= kp;
    }

    out << "SCALARS sound_speed double 1" << nl;
    out << "LOOKUP_TABLE default" << nl;
    p = p1;
    for (label ip = 0; ip < np; ip++)
    {
        T = T1;
        for (label iT = 0; iT < nT; iT++)
        {
            out << gp->c(p,T) << nl;
            T += dT;
        }
        p *= kp;
    }

    out << "SCALARS Cp double 1" << nl;
    out << "LOOKUP_TABLE default" << nl;
    p = p1;
    for (label ip = 0; ip < np; ip++)
    {
        T = T1;
        for (label iT = 0; iT < nT; iT++)
        {
            out << gp->Cp(p,T) << nl;
            T += dT;
        }
        p *= kp;
    }

    out << "SCALARS Cv double 1" << nl;
    out << "LOOKUP_TABLE default" << nl;
    p = p1;
    for (label ip = 0; ip < np; ip++)
    {
        T = T1;
        for (label iT = 0; iT < nT; iT++)
        {
            out << gp->Cv(p,T) << nl;
            T += dT;
        }
        p *= kp;
    }

    out << "SCALARS T_backward double 1" << nl;
    out << "LOOKUP_TABLE default" << nl;
    p = p1;
    for (label ip = 0; ip < np; ip++)
    {
        T = T1;
        for (label iT = 0; iT < nT; iT++)
        {
            scalar hs = gp->Hs(p,T);
            out << gp->THs(hs, p, 300.) << nl;
            T += dT;
        }
        p *= kp;
    }

    
    runTime.printExecutionTime(Info);

    Info << "Created thermo.vtk" << nl;
    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
