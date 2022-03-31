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

#include "subsonicInletTotalFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "totalTemperatureFvPatchScalarField.H"
// #include "inletOutletTotalTemperatureFvPatchScalarField.H"
#include "psiThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::subsonicInletTotalFvPatchVectorField::
subsonicInletTotalFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(p, iF),
    TName_("T"),
    inletDir_(p.size())
{
    refValue() = *this;
    refGrad() = Zero;
    valueFraction() = 0.0;
}


Foam::subsonicInletTotalFvPatchVectorField::
subsonicInletTotalFvPatchVectorField
(
    const subsonicInletTotalFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchVectorField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    inletDir_(mapper(ptf.inletDir_))
{}


Foam::subsonicInletTotalFvPatchVectorField::
subsonicInletTotalFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchVectorField(p, iF),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    inletDir_("inletDirection", dict, p.size())
{
    patchType() = dict.lookupOrDefault<word>("patchType", word::null);
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
    refValue() = *this;
    refGrad() = Zero;
    valueFraction() = 0.0;
}


Foam::subsonicInletTotalFvPatchVectorField::
subsonicInletTotalFvPatchVectorField
(
    const subsonicInletTotalFvPatchVectorField& pivpvf
)
:
    mixedFvPatchVectorField(pivpvf),
    TName_(pivpvf.TName_),
    inletDir_(pivpvf.inletDir_)
{}


Foam::subsonicInletTotalFvPatchVectorField::
subsonicInletTotalFvPatchVectorField
(
    const subsonicInletTotalFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(pivpvf, iF),
    TName_(pivpvf.TName_),
    inletDir_(pivpvf.inletDir_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::subsonicInletTotalFvPatchVectorField::autoMap
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


void Foam::subsonicInletTotalFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    mixedFvPatchVectorField::rmap(ptf, addr);

    const subsonicInletTotalFvPatchVectorField& tiptf =
        refCast<const subsonicInletTotalFvPatchVectorField>
        (ptf);

    inletDir_.rmap(tiptf.inletDir_, addr);
}


void Foam::subsonicInletTotalFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const auto& thermo =
        mesh.lookupObject<psiThermo>("thermophysicalProperties");

    // Need R of the free-stream flow.  Assume R is independent of location
    // along patch so use face 0
    label cellI = patch().faceCells()[0];
    scalar Cp = thermo.Cp()()[cellI];
    scalar Cv = thermo.Cv()()[cellI];
    scalar R = Cp - Cv;
    scalar gamma = Cp / Cv;

    const volScalarField& T = 
        db().lookupObject<volScalarField>(TName_);

    const fvPatchScalarField& pT = 
        patch().patchField<volScalarField, scalar>(T);

    scalarField pT0;
    if (pT.type() == "totalTemperature")
    {
        pT0 = refCast<const totalTemperatureFvPatchScalarField>(pT).T0();
    }
//    else if (pT.type() == "inletOutletTotalTemperature")
//    {
//        pT0 = refCast<const inletOutletTotalTemperatureFvPatchScalarField>(pT).T0();
//    }
    else
    {
        FatalErrorIn("subsonicInletTotalFvPatchVectorField::updateCoeffs()") 
            << "the subsonicInletTotal has to be combined either with "
            << "totalTemperature or inletOutletTotalTemperature condition!"
            << abort(FatalError);
    }

    const vectorField& Uint = internalField();
    const scalarField& Tint = T.internalField();

    const vectorField& pSf = patch().Sf();

    vectorField& refValue = this->refValue();
    scalarField& valFraction = this->valueFraction();

    forAll(pSf, faceI) {
        // Inward normal
        vector n = -pSf[faceI] / mag(pSf[faceI]);
        vector dir = inletDir_[faceI] / mag(inletDir_[faceI]);
	
        label faceCellI = patch().faceCells()[faceI];

        // Total temperature
        scalar T0 = pT0[faceI];
        
        // Normal velocity in the inner cell (inward normal)
        scalar u1 = n & Uint[faceCellI];

        // Sound speed in the inner cell
        scalar c1 = sqrt(gamma*R*Tint[faceCellI]);
	
        // Riemann invariant
        scalar Rm = u1 - 2*c1/(gamma-1);
	
	// Calculate sound speed at the boundary via quadratic eq.
        scalar cb;
        scalar oneByCos = 1 / (n & dir);
        {
	  scalar tmp= oneByCos*oneByCos; 
	  scalar a = 1+2.0/(gamma-1)*tmp;
	  scalar b = 2*tmp*Rm;
	  
	  scalar c = (gamma-1)/2.*tmp*Rm*Rm - gamma*R*T0;
	  //Perr << a << "\t" << b << "\t" << c << "\t" << b*b-4*a*c << endl;
	  cb = (-b+sqrt(max(b*b-4*a*c,0.0)))/(2*a);
        }
        
        scalar ub = 2*cb/(gamma-1) + Rm;
        scalar uMag = ub * oneByCos;
        refValue[faceI] = dir*uMag;
        valFraction[faceI] = pos(ub);
    }
    
    mixedFvPatchVectorField::updateCoeffs();
}


void Foam::subsonicInletTotalFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
#if (OPENFOAM_PLUS>=1712 || OPENFOAM >= 1812)
    os.writeEntryIfDifferent<word>("T", "T", TName_);
 #else
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

void Foam::subsonicInletTotalFvPatchVectorField::operator=
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
        subsonicInletTotalFvPatchVectorField
    );
}

// ************************************************************************* //
