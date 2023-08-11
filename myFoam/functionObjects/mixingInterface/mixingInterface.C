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

\*---------------------------------------------------------------------------*/

#include "mixingInterface.H"
#include "Time.H"
#include "fvCFD.H"
#include "addToRunTimeSelectionTable.H"

#include "totalPressureFvPatchScalarField.H"
#include "totalTemperatureFvPatchScalarField.H"
#include "subsonicInletTotalFvPatchVectorField.H"
#include "myPressureDirectedInletVelocityFvPatchVectorField.H"
#include "isentropicPressureFvPatchScalarField.H"
#include "isentropicTemperatureFvPatchScalarField.H"
#include "isentropicInletVelocityFvPatchVectorField.H"
#include "psiThermo.H"
#include "gasProperties.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(mixingInterface, 0);
    addToRunTimeSelectionTable(functionObject, mixingInterface, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::fvMesh&
Foam::functionObjects::mixingInterface::mesh() const
{
    return refCast<const fvMesh>(obr_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::mixingInterface::mixingInterface
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::regionFunctionObject(name, runTime, dict),
    upstreamPatch_(dict.get<word>("upstreamPatch")),
    downstreamPatch_(dict.get<word>("downstreamPatch")),
    axis_(dict.get<vector>("axis")),
    origin_(dict.getOrDefault<vector>("origin", Zero)),
    configuration_(AXIAL),
    segments_(dict.getOrDefault<label>("segments", 10)),
    convectiveVariables_(dict.getOrDefault<wordList>("convective", wordList())),
    relax_(dict.getOrDefault<scalar>("relaxation", 1.0)),
    frequency_(dict.getOrDefault<label>("frequency", 1)),
    upstreamPatchID_(mesh().boundaryMesh().findPatchID(upstreamPatch_)),
    downstreamPatchID_(mesh().boundaryMesh().findPatchID(downstreamPatch_)),
    counter_(Zero)
{
    //Info << "MIXINGINTERFACE CONSTRUCTOR" << endl;
    read(dict);

    axis_ /= mag(axis_);
    
    word config = dict.get<word>("configuration");
    if (config == "radial")
    {
        configuration_ = RADIAL;
        parameter_ = [&](vector x)
            {
                vector pos = x - origin_;
                return (pos & axis_);
            };
    }
    else if (config == "axial")
    {
        configuration_ = AXIAL;        
        parameter_ = [&](vector x)
            {
                vector pos = x - origin_;
                vector r = pos - (axis_ & pos)*axis_;
                return mag(r);
            };
    }
    else
    {
        FatalError 
            << "Configuration can be either radial or axial!" << nl
            << exit(FatalError);
    }

    
    if (upstreamPatchID_ < 0)
    {
        FatalErrorInFunction
            << "Unable to find upstream patch " << upstreamPatch_ << abort(FatalError);
    }
    if (downstreamPatchID_ < 0)
    {
        FatalErrorInFunction
            << "Unable to find downstream patch " << downstreamPatch_ << abort(FatalError);
    }

    const polyPatch& upPatch   = mesh().boundaryMesh()[upstreamPatchID_];
    const polyPatch& downPatch = mesh().boundaryMesh()[downstreamPatchID_];

    scalar parMin = 1.e32, parMax = -1.e32;
    for (auto x : upPatch.faceCentres())
    {
        scalar par = parameter_(x);
        parMin = min(parMin, par);
        parMax = max(parMax, par);
    }
    for (auto x : downPatch.faceCentres())
    {
        scalar par = parameter_(x);
        parMin = min(parMin, par);
        parMax = max(parMax, par);
    }
    reduce(parMin, minOp<scalar>());
    reduce(parMax, maxOp<scalar>());
    
    paramBins_.resize(segments_ + 1);
    for (label i=0; i<=segments_; i++)
    {
        scalar xi = -M_PI/2 + M_PI*scalar(i)/segments_;
        paramBins_[i] = parMin + (parMax - parMin)*(sin(xi)+1)/2;
        #ifdef FULLDEBUG
          Info << paramBins_[i] << nl;
        #endif
    }

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::mixingInterface::read(const dictionary& dict)
{
    return true;
}


bool Foam::functionObjects::mixingInterface::execute()
{
    counter_++;
    
    if (counter_ % frequency_ != 0) return true;
    
    Info << "MIXING INTERFACE EXECUTE" << nl << nl;;

    const volScalarField& p  = mesh().lookupObject<volScalarField>(word("p"));
    const volScalarField& T  = mesh().lookupObject<volScalarField>(word("T"));
    const volVectorField& U  = mesh().lookupObject<volVectorField>(word("U"));

    const scalarField& pUp = p.boundaryField()[upstreamPatchID_];
    const scalarField& TUp = T.boundaryField()[upstreamPatchID_];
    const vectorField& UUp = U.boundaryField()[upstreamPatchID_];

    const scalarField& pDn = p.boundaryField()[downstreamPatchID_];

    const auto& thermo =
        mesh().lookupObject<psiThermo>("thermophysicalProperties");
    autoPtr<gasProperties> gasProps(gasProperties::New(thermo));
    
    // Total pressure and total temperature to downstream
    {
        scalarField pTot(pUp.size(), Zero);
        scalarField TTot(pUp.size(), Zero);
        
        forAll(pUp, faceI)
        {
            scalar S  = gasProps->S(pUp[faceI], TUp[faceI]);
            scalar H0 = gasProps->Hs(pUp[faceI], TUp[faceI]) + 0.5*magSqr(UUp[faceI]);
            pTot[faceI] = gasProps->pHS(H0, S, pUp[faceI]);
            TTot[faceI] = gasProps->TpS(pTot[faceI], S, TUp[faceI]);
        }
        
        List<scalar> pTotAvgUp = getScalarAverages(
            pTot,
            mesh().magSf().boundaryField()[upstreamPatchID_],
            mesh().Cf().boundaryField()[upstreamPatchID_]
        );

        List<scalar> TTotAvgUp = getScalarAverages(
            TTot,
            mesh().magSf().boundaryField()[upstreamPatchID_],
            mesh().Cf().boundaryField()[upstreamPatchID_]
        );

        volScalarField::Boundary& pb =
            const_cast<volScalarField*>(&p)->boundaryFieldRef();
        scalarField* p0Ptr = NULL;

        if (pb[downstreamPatchID_].type() == "totalPressure")
        {
            auto& pp = refCast<totalPressureFvPatchScalarField>(pb[downstreamPatchID_]);
            p0Ptr = &pp.p0();
        }
        else if (pb[downstreamPatchID_].type() == "isentropicPressure")
        {
            auto& pp = refCast<isentropicPressureFvPatchScalarField>(pb[downstreamPatchID_]);
            p0Ptr = &pp.p0();
        }
	else 
	{
            FatalErrorIn("mixingInterface::execute()") 
                << "the mixing interface works only with totalPressure "
                << "or isentropicPressure boundary conditions!"
                << abort(FatalError);
        }        
        auto& p0 = *p0Ptr;

        forAll(p0, i)
        {
            vector x = mesh().Cf().boundaryField()[downstreamPatchID_][i];
            p0[i] = (1 - relax_)*p0[i] + relax_*interpolateScalar(pTotAvgUp, parameter_(x));
        }

        volScalarField::Boundary& Tb =
            const_cast<volScalarField*>(&T)->boundaryFieldRef();

        scalarField* T0Ptr = NULL;

        if (Tb[downstreamPatchID_].type() == "totalTemperature")
        {
            auto& pT = refCast<totalTemperatureFvPatchScalarField>(Tb[downstreamPatchID_]);
            T0Ptr = &pT.T0();
        }
        else if (Tb[downstreamPatchID_].type() == "isentropicTemperature")
        {
            auto& pT = refCast<isentropicTemperatureFvPatchScalarField>(Tb[downstreamPatchID_]);
            T0Ptr = &pT.T0();
        }
	else 
	{
            FatalErrorIn("mixingInterface::execute()") 
                << "the mixing interface works only with totalTemperature "
                << "or isentropicTemperature boundary conditions!"
                << abort(FatalError);
        }
        
        auto& T0 = *T0Ptr;
        forAll(T0, i)
        {
            vector x = mesh().Cf().boundaryField()[downstreamPatchID_][i];
            T0[i] = (1 - relax_)*T0[i] + relax_*interpolateScalar(TTotAvgUp, parameter_(x));
        }
    }
        
    // Direction of the velocity downstream
    {
        List<vector> UAvgUp = getVectorAverages(
            UUp,
            mesh().magSf().boundaryField()[upstreamPatchID_],
            mesh().Cf().boundaryField()[upstreamPatchID_]
        );

        volVectorField::Boundary& Ub =
            const_cast<volVectorField*>(&U)->boundaryFieldRef();

	vectorField* inletDirPtr = NULL;
	if (Ub[downstreamPatchID_].type() == "subsonicInletTotal")
	{
            auto& pU = refCast<subsonicInletTotalFvPatchVectorField>(Ub[downstreamPatchID_]);
	    inletDirPtr = &pU.inletDirection();
	}
	else if (Ub[downstreamPatchID_].type() == "myPressureDirectedInletVelocity")
	{
            auto& pU = refCast<myPressureDirectedInletVelocityFvPatchVectorField>(Ub[downstreamPatchID_]);
	    inletDirPtr = &pU.inletDirection();
	}
	else if (Ub[downstreamPatchID_].type() == "isentropicInletVelocity")
	{
            auto& pU = refCast<isentropicInletVelocityFvPatchVectorField>(Ub[downstreamPatchID_]);
	    inletDirPtr = &pU.inletDirection();
	}
	else 
	{
            FatalErrorIn("mixingInterface::execute()") 
                << "the mixing interface works only with subsonicInletTotal, "
                << "isentropicInletVelocity or myPressureDirectedInletVelocity "
                << "boundary conditions!"
                << abort(FatalError);
        }
      	auto& inletDir = *inletDirPtr;

        forAll(inletDir, i)
        {
            vector x = mesh().Cf().boundaryField()[downstreamPatchID_][i];
            vector v = interpolateVector(UAvgUp, parameter_(x));

            vector r = (x - origin_) - (axis_ & (x - origin_))*axis_;
            r /= mag(r);
            vector t = axis_ ^ r;
            
            vector vxy = v[0]*axis_ + v[1]*r + v[2]*t;
            vxy /= max(mag(vxy), 1.e-10);
            
            inletDir[i] = (1 - relax_)*inletDir[i] + relax_*vxy;
            inletDir[i] /= mag(inletDir[i]);

            // Check if the direction goes into the domain
            vector n = mesh().Sf().boundaryField()[downstreamPatchID_][i];
            n /= mag(n);
            if ( (inletDir[i] & n) >= 0)
            {
                inletDir[i] -= 2*(inletDir[i] & n)*n;
            }
        }
    }


    // Static pressure to upstream (using mean values)
    {
        List<scalar> pAvgDn = getScalarAverages(
            pDn,
            mesh().magSf().boundaryField()[downstreamPatchID_],
            mesh().Cf().boundaryField()[downstreamPatchID_]
        );

        volScalarField::Boundary& pb =
            const_cast<volScalarField*>(&p)->boundaryFieldRef();

        fixedValueFvPatchScalarField& pp =
            refCast<fixedValueFvPatchScalarField>(pb[upstreamPatchID_]);

        // Extrapolate intrnal pressure 
        pp.operator=(pp.patchInternalField());
        
        List<scalar> pAvgUp = getScalarAverages(
            pUp,
            mesh().magSf().boundaryField()[upstreamPatchID_],
            mesh().Cf().boundaryField()[upstreamPatchID_]
        );

        List<scalar> dp = pAvgDn - pAvgUp;
        
        forAll(pp, i)
        {
            vector x = mesh().Cf().boundaryField()[upstreamPatchID_][i];
            //pp[i] = (1 - relax_)*pp[i] + relax_*interpolateScalar(pAvgDn, parameter_(x));
            pp[i] += relax_*interpolateScalar(dp, parameter_(x));
        }
    }


    // Convective fields to downstream
    for (auto name : convectiveVariables_)
    {
        const volScalarField& v  = mesh().lookupObject<volScalarField>(name);
        const fvPatchScalarField& vUp = v.boundaryField()[upstreamPatchID_];

        List<scalar> vAvgDn = getScalarAverages(
            vUp,
            mesh().magSf().boundaryField()[upstreamPatchID_],
            mesh().Cf().boundaryField()[upstreamPatchID_]
        );

        volScalarField::Boundary& vb =
            const_cast<volScalarField*>(&v)->boundaryFieldRef();

        fixedValueFvPatchScalarField& pv =
            refCast<fixedValueFvPatchScalarField>(vb[downstreamPatchID_]);
        forAll(pv, i)
        {
            vector x = mesh().Cf().boundaryField()[downstreamPatchID_][i];
            pv[i] = (1 - relax_)*pv[i] + relax_*interpolateScalar(vAvgDn, parameter_(x));
        }        
    }
    
    return true;
}


bool Foam::functionObjects::mixingInterface::end()
{
    return true;
}


bool Foam::functionObjects::mixingInterface::write()
{
    return true;
}


Foam::List<scalar>  Foam::functionObjects::mixingInterface::getScalarAverages(
    const scalarField& values,
    const scalarField& weights,
    const vectorField& coords
) const
{
    List<scalar> averages(segments_);
    List<scalar> w(segments_);

    for (label i=0; i<segments_; i++)
    {
        averages[i] = 0.0;
        w[i] = 0.0;
    };

    forAll(coords, i)
    {
        scalar par = parameter_(coords[i]);
        
        label bin = getBin(par);
        averages[bin] += weights[i]*values[i];
        w[bin] += weights[i];
    }

    reduce(averages, sumOp<scalarList>());
    reduce(w, sumOp<scalarList>());

    for (label i=0; i<segments_; i++)
    {
        if (w[i]<1.e-16)
            Info << "ZERO WEIGHTS!" << nl;
        averages[i] /= w[i];
    };

    return averages;
}


Foam::List<vector>  Foam::functionObjects::mixingInterface::getVectorAverages(
    const vectorField& values,
    const scalarField& weights,
    const vectorField& coords
) const
{
    scalarField axf = values & axis_;
    List<scalar> axialAverages = getScalarAverages(axf, weights, coords);

    vectorField r = (coords - origin_) - (axis_ & (coords - origin_))*axis_;
    r /= mag(r);
    
    scalarField rf = values & r;
    List<scalar> radialAverages = getScalarAverages(rf, weights, coords);

    scalarField tf = values & (axis_ ^ r);
    List<scalar> tangentialAverages = getScalarAverages(tf, weights, coords);
    
    List<vector> averages(segments_);
    for (label i=0; i<segments_; i++)
    {
        averages[i] = vector(axialAverages[i], radialAverages[i], tangentialAverages[i]);
    }
    
    return averages;
}

Foam::label Foam::functionObjects::mixingInterface::getBin(scalar p) const
{
    label i;
    for (i=1; i<=segments_; i++)
        if (p < paramBins_[i]) return i-1;
    return segments_-1;
}

Foam::scalar Foam::functionObjects::mixingInterface::interpolateScalar(
    const List<scalar>& y,
    scalar p
) const
{
    return y[getBin(p)];
}

Foam::vector Foam::functionObjects::mixingInterface::interpolateVector(
    const List<vector>& y,
    scalar p
) const
{
    return y[getBin(p)];
}


// ************************************************************************* //
