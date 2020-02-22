/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "meanTotalPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meanTotalPressureFvPatchScalarField::meanTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_("U"),
    phiName_("phi"),
    rhoName_("none"),
    psiName_("none"),
    gamma_(0.0),
    p0_(p.size(), 0.0),
    meanValue_(0.0)
{}


Foam::meanTotalPressureFvPatchScalarField::meanTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "none")),
    psiName_(dict.lookupOrDefault<word>("psi", "none")),
    gamma_(psiName_ != "none" ? readScalar(dict.lookup("gamma")) : 1),
    p0_("p0", dict, p.size()),
    meanValue_(readScalar(dict.lookup("meanValue")))
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(p0_);
    }
}


Foam::meanTotalPressureFvPatchScalarField::meanTotalPressureFvPatchScalarField
(
    const meanTotalPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    psiName_(ptf.psiName_),
    gamma_(ptf.gamma_),
    p0_(mapper(ptf.p0_)),
    meanValue_(ptf.meanValue_)
{}


Foam::meanTotalPressureFvPatchScalarField::meanTotalPressureFvPatchScalarField
(
    const meanTotalPressureFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    psiName_(tppsf.psiName_),
    gamma_(tppsf.gamma_),
    p0_(tppsf.p0_),
    meanValue_(tppsf.meanValue_)
{}


Foam::meanTotalPressureFvPatchScalarField::meanTotalPressureFvPatchScalarField
(
    const meanTotalPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    psiName_(tppsf.psiName_),
    gamma_(tppsf.gamma_),
    p0_(tppsf.p0_),
    meanValue_(tppsf.meanValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meanTotalPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    m(p0_, p0_);
}


void Foam::meanTotalPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const meanTotalPressureFvPatchScalarField& tiptf =
        refCast<const meanTotalPressureFvPatchScalarField>(ptf);

    p0_.rmap(tiptf.p0_, addr);
}


void Foam::meanTotalPressureFvPatchScalarField::updateCoeffs
(
    scalarField& p0p,
    const vectorField& Up
)
{
    if (updated())
    {
        return;
    }


    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    if (internalField().dimensions() == dimPressure)
    {
        if (psiName_ == "none")
        {
            // Variable density and low-speed compressible flow
            
            const fvPatchField<scalar>& rho =
                patch().lookupPatchField<volScalarField, scalar>(rhoName_);
            
            operator==(p0p - 0.5*rho*(1.0 - pos(phip))*magSqr(Up));
        }
        else
        {
            // High-speed compressible flow

            const fvPatchField<scalar>& psip =
                patch().lookupPatchField<volScalarField, scalar>(psiName_);
            
            if (gamma_ > 1)
            {
                scalar gM1ByG = (gamma_ - 1)/gamma_;

                scalar relax = 0.1;
                
                p0p = (1-relax) * p0p 
                    + relax * this->patchInternalField()
                    * pow
                    (
                        (1.0 + 0.5*psip*gM1ByG*(1.0 - pos(phip))*magSqr(Up)),
                        1.0/gM1ByG
                    );
                
                scalar pTotAvg = gSum(this->patch().magSf() * p0p) /
                    gSum(this->patch().magSf());
                
                if (meanValue_ < pTotAvg )
                    p0p *= meanValue_ / pTotAvg;
                else
                    p0p += (meanValue_ - pTotAvg);
                
                Info << "pTotAvg = " << pTotAvg << endl;
                
                operator==
                (
                    p0p
                   /pow
                    (
                        (1.0 + 0.5*psip*gM1ByG*(1.0 - pos(phip))*magSqr(Up)),
                        1.0/gM1ByG
                    )
                );
            }
            else
            {
                operator==(p0p/(1.0 + 0.5*psip*(1.0 - pos(phip))*magSqr(Up)));
            }
        }
        
    }
    else if (internalField().dimensions() == dimPressure/dimDensity)
    {
        // Incompressible flow
        operator==(p0p - 0.5*(1.0 - pos(phip))*magSqr(Up));
    }
    
    else
    {
        FatalErrorIn
        (
            "meanTotalPressureFvPatchScalarField::updateCoeffs()"
        )   << " rho or psi set inconsistently, rho = " << rhoName_
            << ", psi = " << psiName_ << ".\n"
            << "    Set either rho or psi or neither depending on the "
               "definition of total pressure." << nl
            << "    Set the unused variable(s) to 'none'.\n"
            << "    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::meanTotalPressureFvPatchScalarField::updateCoeffs()
{
    updateCoeffs
    (
        p0(),
        patch().lookupPatchField<volVectorField, vector>(UName())
    );
}


void Foam::meanTotalPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    #if OPENFOAM_PLUS>=1712
    os.writeEntryIfDifferent<word>("U", "U", UName_);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    #else
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    #endif
    os.writeKeyword("rho") << rhoName_ << token::END_STATEMENT << nl;
    os.writeKeyword("psi") << psiName_ << token::END_STATEMENT << nl;
    os.writeKeyword("gamma") << gamma_ << token::END_STATEMENT << nl;
    writeEntry(os, "p0", p0_);
    os.writeKeyword("meanValue") << meanValue_ << token::END_STATEMENT << nl;
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        meanTotalPressureFvPatchScalarField
    );
}

// ************************************************************************* //
