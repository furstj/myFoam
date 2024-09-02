/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "isentropicInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "isentropicTemperatureFvPatchScalarField.H"
#include "isentropicPressureFvPatchScalarField.H"
#include "psiThermo.H"
#include "gasProperties.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isentropicInletVelocityFvPatchVectorField::
isentropicInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(p, iF),
    pName_("p"),
    TName_("T")
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = Zero;
}


Foam::isentropicInletVelocityFvPatchVectorField::
isentropicInletVelocityFvPatchVectorField
(
    const isentropicInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchVectorField(ptf, p, iF, mapper),
    pName_(ptf.pName_),
    TName_(ptf.TName_)
{
    if (ptf.inletDir_.size())
    {
        inletDir_ = mapper(ptf.inletDir_);
    }
    if (ptf.tangentialVelocity_.size())
    {
        tangentialVelocity_ = mapper(ptf.tangentialVelocity_);
    }
}


Foam::isentropicInletVelocityFvPatchVectorField::
isentropicInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchVectorField(p, iF),
    pName_(dict.getOrDefault<word>("p", "p")),
    TName_(dict.getOrDefault<word>("T", "T"))
{
    bool hasTangentialVelocity = false;
    patchType() = dict.getOrDefault<word>("patchType", word::null);
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
    if (dict.found("tangentialVelocity"))
    {
        setTangentialVelocity
        (
            vectorField("tangentialVelocity", dict, p.size())
        );
        hasTangentialVelocity = true;
    }
    if (dict.found("inletDirection"))
    {
        if (hasTangentialVelocity)
        {
            FatalIOErrorInFunction(dict)
                << "For " << this->internalField().name() << " on "
                << this->patch().name() << nl
                << "Doesn't allow both: inletDirection and tangentialVelocity" << nl
                << exit(FatalIOError);
        }
        
    }
    refValue() = *this;
    refGrad() = Zero;
    valueFraction() = Zero;
}


Foam::isentropicInletVelocityFvPatchVectorField::
isentropicInletVelocityFvPatchVectorField
(
    const isentropicInletVelocityFvPatchVectorField& pivpvf
)
:
    mixedFvPatchVectorField(pivpvf),
    pName_(pivpvf.pName_),
    TName_(pivpvf.TName_),
    inletDir_(pivpvf.inletDir_),
    tangentialVelocity_(pivpvf.tangentialVelocity_)
{}


Foam::isentropicInletVelocityFvPatchVectorField::
isentropicInletVelocityFvPatchVectorField
(
    const isentropicInletVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(pivpvf, iF),
    pName_(pivpvf.pName_),
    TName_(pivpvf.TName_),
    inletDir_(pivpvf.inletDir_),
    tangentialVelocity_(pivpvf.tangentialVelocity_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::isentropicInletVelocityFvPatchVectorField::
setTangentialVelocity(const vectorField& tangentialVelocity)
{
    tangentialVelocity_ = tangentialVelocity;
    const vectorField n(patch().nf());
    refValue() = tangentialVelocity_ - n*(n & tangentialVelocity_);
}

void Foam::isentropicInletVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchVectorField::autoMap(m);
#if (OPENFOAM >= 1812)
    if (inletDir_.size())
    {
        inletDir_.autoMap(m);
    }
    if (tangentialVelocity_.size())
    {
        tangentialVelocity_.autoMap(m);
    }
#else
    if (inletDir_.size())
    {
        m(inletDir_, inletDir_);
    }
    if (tangentialVelocity_.size())
    {
        m(tangentialVelocity_, tangentialVelocity_);
    }
#endif    
}


void Foam::isentropicInletVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    mixedFvPatchVectorField::rmap(ptf, addr);

    const isentropicInletVelocityFvPatchVectorField& tiptf =
        refCast<const isentropicInletVelocityFvPatchVectorField>
        (ptf);
    if (inletDir_.size())
    {
        inletDir_.rmap(tiptf.inletDir_, addr);
    }
    if (tangentialVelocity_.size())
    {
        tangentialVelocity_.rmap(tiptf.tangentialVelocity_, addr);
    }
    
}


void Foam::isentropicInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const auto& thermo =
        mesh.lookupObject<psiThermo>("thermophysicalProperties");

    autoPtr<gasProperties> gasProps(gasProperties::New(thermo));
    
    const volScalarField& T = 
        db().lookupObject<volScalarField>(TName_);

    const fvPatchScalarField& pT = 
        patch().patchField<volScalarField, scalar>(T);

    const volScalarField& p = 
        db().lookupObject<volScalarField>(pName_);

    const fvPatchScalarField& pp = 
        patch().patchField<volScalarField, scalar>(p);
    
    scalarField pT0, pp0;
    if (pT.type() == "isentropicTemperature" && pp.type() == "isentropicPressure")
    {
        pT0 = refCast<const isentropicTemperatureFvPatchScalarField>(pT).T0();
        pp0 = refCast<const isentropicPressureFvPatchScalarField>(pp).p0();
    }
    else
    {
        FatalErrorIn("isentropicInletVelocityFvPatchVectorField::updateCoeffs()") 
            << "the isentropicInletVelocity has to be combined with "
            << "isentropicTemperature and isentropicPressure conditions!"
            << abort(FatalError);
    }

    const scalarField& pint = p.internalField();
    const scalarField& Tint = T.internalField();
    const vectorField& Uint = internalField();
    
    const vectorField& pSf = patch().Sf();

    vectorField& refValue = this->refValue();
    scalarField& valFraction = this->valueFraction();

    const bool hasInletDirection = inletDir_.size() > 0;
    const bool hasTangentialVelocity = tangentialVelocity_.size() > 0;

    forAll(pSf, faceI) {	
        label faceCellI = patch().faceCells()[faceI];

        // Total temperature, total pressure, ...
        const scalar T0 = pT0[faceI];
        const scalar p0 = pp0[faceI];
        const scalar S  = gasProps->S(p0, T0);
        const scalar H0 = gasProps->Hs(p0, T0);

        /*
        // Pressure extrapolation
        const scalar pb = pint[faceCellI];
        const scalar Tb = gasProps->TpS(pb, S, pT[faceI]);
        const scalar hb = gasProps->Hs(pb, Tb);
        const scalar magU = sqrt(2*max(H0 - hb, 0.0));
        */

        // Inward normal
        const vector n = -pSf[faceI] / mag(pSf[faceI]);
        const vector dir = (hasInletDirection) ? inletDir_[faceI] / mag(inletDir_[faceI]) : n;
        const vector utau = (hasTangentialVelocity) ? tangentialVelocity_[faceI] : Zero;

        const scalar oneByCos = 1 / (n & dir);
        const scalar p1 = pint[faceCellI];
        const scalar u1 = n & Uint[faceCellI];
        const scalar c1 = gasProps->c(pint[faceCellI], Tint[faceCellI]);
        const scalar rho1 = gasProps->rho(pint[faceCellI], Tint[faceCellI]);
        
        scalar pb = pp[faceI];
        scalar Tb = pT[faceI];
        label iter = 0;
        const label maxIter = 100;
        const scalar pTol = 1.e-6*pb;
        scalar dp, dH, ub;
        do
        {
            ub = u1 + (pb - p1)/(rho1*c1);
            ub = max(ub, 0.0);
            Tb = gasProps->TpS(pb, S, Tb);
            dH = H0 - gasProps->Hs(pb, Tb) - 0.5*(sqr(ub*oneByCos) + magSqr(utau));
            dp = dH/(1/gasProps->rho(pb,Tb) + ub*sqr(oneByCos)/(rho1*c1));
            pb += dp;
            //Info << iter << " ub  " << ub << "  dp  " << dp << nl;
            if (iter++ > maxIter)
            {
                FatalErrorInFunction
                    << "Maximum number of iterations exceeded: " << maxIter
                        << " when starting from p:" << pp[faceI] << nl
                        << " T  : " << Tb << nl
                        << " p  : " << pb << nl
                        << " ub  : " << ub << nl
                        << " dp  : " << dp << nl
                        << " dH  : " << dH << nl
                        << " tol: " << pTol
                        << abort(FatalError);
            }
        } while (mag(dp) > pTol);
        
        refValue[faceI] = dir*ub*oneByCos + utau;
        valFraction[faceI] = pos0(ub);
    }
    
    mixedFvPatchVectorField::updateCoeffs();
}


void Foam::isentropicInletVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
#if (OPENFOAM_PLUS>=1712 || OPENFOAM >= 1812)
    os.writeEntryIfDifferent<word>("p", "p", pName_);
    os.writeEntryIfDifferent<word>("T", "T", TName_);
 #else
    writeEntryIfDifferent<word>(os, "p", "p", pName_);
    writeEntryIfDifferent<word>(os, "T", "T", TName_);
 #endif
#if (OPENFOAM >= 1812)
    if (inletDir_.size())
    {
        inletDir_.writeEntry("inletDirection", os);
    }
    if (tangentialVelocity_.size())
    {
        tangentialVelocity_.writeEntry("tangentialVelocity", os);
    }
    this->writeEntry("value", os);
#else
    if (inletDir_.size())
    {
        writeEntry(os, "inletDirection", inletDir_);
    }
    if (tangentialVelocity_.size())
    {
        writeEntry(os, "tangentialVelocity", tangentialVelocity_);
    }
    
    writeEntry(os, "value", *this);
#endif    

}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        isentropicInletVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
