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
    TName_("T"),
    inletDir_(p.size())
{
    refValue() = *this;
    refGrad() = Zero;
    valueFraction() = 0.0;
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
    TName_(ptf.TName_),
    inletDir_(ptf.inletDir_, mapper)
{}


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
    TName_(dict.getOrDefault<word>("T", "T")),
    inletDir_("inletDirection", dict, p.size())
{
    patchType() = dict.getOrDefault<word>("patchType", word::null);
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
    refValue() = *this;
    refGrad() = Zero;
    valueFraction() = 0.0;
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
    inletDir_(pivpvf.inletDir_)
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
    inletDir_(pivpvf.inletDir_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::isentropicInletVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchVectorField::autoMap(m);
#if (OPENFOAM >= 1812)
    inletDir_.autoMap(m);
#else
    m(inletDir_, inletDir_);
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

    inletDir_.rmap(tiptf.inletDir_, addr);
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

    const vectorField& Uint = internalField();
    const scalarField& Tint = T.internalField();
    const scalarField& pint = p.internalField();

    const vectorField& pSf = patch().Sf();

    vectorField& refValue = this->refValue();
    scalarField& valFraction = this->valueFraction();

    forAll(pSf, faceI) {
        // Inward normal
        const vector n = -pSf[faceI] / mag(pSf[faceI]);
        const vector dir = inletDir_[faceI] / mag(inletDir_[faceI]);
	
        label faceCellI = patch().faceCells()[faceI];

        // Total temperature, total pressure, ...
        const scalar T0 = pT0[faceI];
        const scalar p0 = pp0[faceI];
        const scalar S  = gasProps->S(p0, T0);
        const scalar H0 = gasProps->Hs(p0, T0);
        
        // Normal velocity in the inner cell (inward normal)
        const scalar u1 = n & Uint[faceCellI];

        // Sound speed and density in the inner cell
        const scalar p1 = pint[faceCellI];
        const scalar T1 = Tint[faceCellI]; 
        const scalar c1   = gasProps->c(p1, T1);
        const scalar rho1 = gasProps->rho(p1, T1);

        // Calculation using Riemann invariant:
        // 0 = dRminus = (ub - u1) - \int_{\rho_1}^{\rho_b} c(\rho,S)/rho \, d\rho
        //
        // integral is approximated using trapezoidal rule:
        // u_b = u_1 - 0.5*(c1/rho1 + cb/rhob)*(rho1 - rhob)

        scalar pb = pp[faceI];
        scalar Tb = pT[faceI];
        scalar cb = gasProps->c(pb, Tb);
        scalar rhob = gasProps->rho(pb, Tb);
        scalar oneByCos = 1 / (n & dir);
        
        scalar ub = u1;
        scalar dp;

        const label maxIter = 100;
        label iter = 0;
        scalar pTol = 1.e-5*p1;

        scalarList dphist(maxIter);
        
        do
        {
            ub = u1 - 0.5*(c1/rho1 + cb/rhob)*(rho1 - rhob);
            scalar h = H0 - 0.5*sqr(ub*oneByCos);
            scalar dh = h - gasProps->Hs(pb, Tb);
            dp = 0.5*rhob*dh;
            pb += dp;
            Tb = gasProps->TpS(pb, S, Tb);
            rhob = gasProps->rho(pb, Tb);
            cb = gasProps->c(pb, Tb);
            dphist[iter] = dp;
            if (iter++ > maxIter)
            {
                FatalErrorInFunction
                    << "Maximum number of iterations exceeded: " << maxIter << nl
                        << " T  : " << T1 << "  ... " << Tb << nl
                        << " p  : " << p1 << "  ... " << pb << nl
                        << " ub : " << u1 << "  ... " << ub << nl
                        << " S  : " << S << nl
                        << " h  : " << h << "  dh   " << dh << nl
                        << " p  : " << pb << "  dp   " << dp << nl
                        << " history  : " << dphist << nl
                        << abort(FatalError);
            }           
        } while (mag(dp) > pTol);
        
        scalar uMag = ub * oneByCos;
        refValue[faceI] = dir*uMag;
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
    inletDir_.writeEntry("inletDirection", os);
    this->writeEntry("value", os);
#else
    writeEntry(os, "inletDirection", inletDir_);
    writeEntry(os, "value", *this);
#endif    

}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::isentropicInletVelocityFvPatchVectorField::operator=
(
    const fvPatchField<vector>& pvf
)
{
    fvPatchField<vector>::operator=
    (
        valueFraction()*(inletDir_*(inletDir_ & pvf))
      + (1 - valueFraction())*pvf
    );
}


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
