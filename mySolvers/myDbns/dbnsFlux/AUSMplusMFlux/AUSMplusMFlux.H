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

Application
    ausmPlusMFlux

Description
    AUSM+M flux scheme. From:

    "An improved AUSM-family scheme with robustness and accuracy for all Mach number flows",
    Chen, Shu-sheng and Cai, Fang-jie and Xue, Hai-chao and  Wang, Ning and Yan, Chao,
    Applied Mathematical Modelling (77), doi.org/10.1016/j.apm.2019.09.005, 2020.

Author(s)
    Adam Tater         - Major development.
    Maximilian Maigler - Minor improvements and testing.

SourceFiles
    AUSMplusMFlux.H

\*---------------------------------------------------------------------------*/

#ifndef AUSMplusMFlux_H
#define AUSMplusMFlux_H

#include "dbnsFlux.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class AUSMplusMFlux Declaration
\*---------------------------------------------------------------------------*/

class AUSMplusMFlux: public dbnsFlux
{
private:
    AUSMplusMFlux();
    AUSMplusMFlux(const AUSMplusMFlux&);

    // Boolean flag for pressure sensing
    bool requiresPressureSensing_;
  
    scalar sqrMachInf_;
    scalar alpha0_;

public:
    TypeName("AUSMplusMFlux");

    // Constructors

    //- Constructor 
    AUSMplusMFlux(const fvMesh&, const dictionary&);


    // Destructor - default
    virtual ~AUSMplusMFlux() {};

    // Member Functions
    
    //- Update fluxes based on current left and right states
    virtual void evaluateFlux
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
    ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
