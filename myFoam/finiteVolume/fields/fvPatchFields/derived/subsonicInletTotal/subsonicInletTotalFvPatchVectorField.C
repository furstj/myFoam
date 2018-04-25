/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
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

#include "subsonicInletTotalFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "totalPressureFvPatchScalarField.H"
#include "totalTemperatureFvPatchScalarField.H"
#include <stdlib.h>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::subsonicInletTotalFvPatchVectorField::
subsonicInletTotalFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    pName_("p"),
    TName_("T"),
    inletDir_(p.size())
{}


Foam::subsonicInletTotalFvPatchVectorField::
subsonicInletTotalFvPatchVectorField
(
    const subsonicInletTotalFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    pName_(ptf.pName_),
    TName_(ptf.TName_),
    inletDir_(ptf.inletDir_, mapper)
{}


Foam::subsonicInletTotalFvPatchVectorField::
subsonicInletTotalFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    pName_(dict.lookupOrDefault<word>("p", "p")),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    inletDir_("inletDirection", dict, p.size())
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}

Foam::subsonicInletTotalFvPatchVectorField::
subsonicInletTotalFvPatchVectorField
(
    const subsonicInletTotalFvPatchVectorField& sfspvf
)
:
    fixedValueFvPatchVectorField(sfspvf),
    pName_(sfspvf.pName_),
    TName_(sfspvf.TName_),
    inletDir_(sfspvf.inletDir_)
{}


Foam::subsonicInletTotalFvPatchVectorField::
subsonicInletTotalFvPatchVectorField
(
    const subsonicInletTotalFvPatchVectorField& sfspvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(sfspvf, iF),
    pName_(sfspvf.pName_),
    TName_(sfspvf.TName_),
    inletDir_(sfspvf.inletDir_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::subsonicInletTotalFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
    inletDir_.autoMap(m);
}


void Foam::subsonicInletTotalFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const subsonicInletTotalFvPatchVectorField& tiptf =
        refCast<const subsonicInletTotalFvPatchVectorField>(ptf);

    inletDir_.rmap(tiptf.inletDir_, addr);
}

void Foam::subsonicInletTotalFvPatchVectorField::updateCoeffs()
{
    if (!size() || updated())
    {
        return;
    }

    const totalTemperatureFvPatchScalarField& pT =
      refCast<const totalTemperatureFvPatchScalarField>
      ( patch().lookupPatchField<volScalarField, scalar>( TName() ) );

    const totalPressureFvPatchScalarField& pp =
      refCast<const totalPressureFvPatchScalarField>
      ( patch().lookupPatchField<volScalarField, scalar>( pName() ) );
							

    const fvPatchField<scalar>& ppsi =
        patch().lookupPatchField<volScalarField, scalar>("thermo:psi");

    // Need R of the free-stream flow.  Assume R is independent of location
    // along patch so use face 0
    scalar R = 1.0/(ppsi[0]*pT[0]);

    scalar gamma = 
        static_cast<const totalPressureFvPatchScalarField*>(&pp)->gamma();

    const scalarField& pT0 = 
        static_cast<const totalTemperatureFvPatchScalarField*>(&pT)->T0();


    const Field<vector>& Uint = internalField();
    const Field<scalar>& Tint =
                   db().lookupObject<volScalarField>( TName() );

    const vectorField& pSf = patch().Sf();

    vectorField& Up = *this;

    forAll(pSf, faceI) {
        // Outward normal
        vector n = -pSf[faceI] / mag(pSf[faceI]);
        vector dir = inletDir()[faceI] / mag(inletDir()[faceI]);

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
            //Info << a << "\t" << b << "\t" << c << endl;
            cb = (-b+sqrt(b*b-4*a*c))/(2*a);
        }
        
        scalar ub = 2*cb/(gamma-1) + Rm;
        scalar uMag = ub * oneByCos;
        Up[faceI] = dir*uMag;

   }
    //std::exit(1);
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::subsonicInletTotalFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    #if OPENFOAM_PLUS>=1712
    os.writeEntryIfDifferent<word>("p", "p", pName_);
    os.writeEntryIfDifferent<word>("T", "T", TName_);
    #else
    writeEntryIfDifferent<word>(os, "p", "p", pName_);
    writeEntryIfDifferent<word>(os, "T", "T", TName_);
    #endif
    inletDir_.writeEntry("inletDirection", os);
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::subsonicInletTotalFvPatchVectorField::operator=
(
    const fvPatchField<vector>& pvf
)
{
    fvPatchField<vector>::operator=(inletDir_*(inletDir_ & pvf));
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
