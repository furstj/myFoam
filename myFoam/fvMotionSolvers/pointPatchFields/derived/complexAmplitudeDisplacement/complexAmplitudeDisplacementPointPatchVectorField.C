/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original authors
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

#include "complexAmplitudeDisplacementPointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "fvMesh.H"
#include "IFstream.H"
#include "transformField.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

complexAmplitudeDisplacementPointPatchVectorField::
complexAmplitudeDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    frequency_(0.0),
    phase_(0.0),
    ramp_(0.0),
    realAmplitude_(p.size(), vector(0,0,0)),
    imagAmplitude_(p.size(), vector(0,0,0))
{}


complexAmplitudeDisplacementPointPatchVectorField::
complexAmplitudeDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    frequency_(readScalar(dict.lookup("frequency"))),
    phase_(readScalar(dict.lookup("phase"))),
    ramp_(readScalar(dict.lookup("ramp"))),
    realAmplitude_("realAmplitude", dict, p.size()),
    imagAmplitude_("imagAmplitude", dict, p.size())
{
    updateCoeffs();
}


complexAmplitudeDisplacementPointPatchVectorField::
complexAmplitudeDisplacementPointPatchVectorField
(
    const complexAmplitudeDisplacementPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    frequency_(ptf.frequency_),
    phase_(ptf.phase_),
    ramp_(ptf.ramp_),
    realAmplitude_(mapper(ptf.realAmplitude_)),
    imagAmplitude_(mapper(ptf.imagAmplitude_))
{}


complexAmplitudeDisplacementPointPatchVectorField::
complexAmplitudeDisplacementPointPatchVectorField
(
    const complexAmplitudeDisplacementPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    frequency_(ptf.frequency_),
    phase_(ptf.phase_),
    ramp_(ptf.ramp_),
    realAmplitude_(ptf.realAmplitude_),
    imagAmplitude_(ptf.imagAmplitude_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void complexAmplitudeDisplacementPointPatchVectorField::autoMap
(
    const pointPatchFieldMapper& m
)
{
    fixedValuePointPatchField<vector>::autoMap(m);
#if (OPENFOAM >= 1812)
    realAmplitude_.autoMap(m);
    imagAmplitude_.autoMap(m);
#else
    m(realAmplitude_, realAmplitude_);
    m(imagAmplitude_, imagAmplitude_);
#endif
}


void complexAmplitudeDisplacementPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    const complexAmplitudeDisplacementPointPatchVectorField& cAptf =
        refCast<const complexAmplitudeDisplacementPointPatchVectorField>(ptf);

    fixedValuePointPatchField<vector>::rmap(cAptf, addr);
    realAmplitude_.rmap(cAptf.realAmplitude_, addr);
    imagAmplitude_.rmap(cAptf.imagAmplitude_, addr);
}


void complexAmplitudeDisplacementPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = this->internalField().mesh()();
    const Time& t = mesh.time();

    scalar phi(2*M_PI*frequency_ * t.value() + phase_);
    scalar K(1.0);
    if (t.value() < ramp_) K = t.value() / ramp_;

    Field<vector>::operator=(K*(realAmplitude_*cos(phi) - imagAmplitude_*sin(phi)));

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void complexAmplitudeDisplacementPointPatchVectorField::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);
    os.writeKeyword("frequency") << frequency_ << token::END_STATEMENT << nl;
    os.writeKeyword("phase") << phase_ << token::END_STATEMENT << nl;
    os.writeKeyword("ramp") << ramp_ << token::END_STATEMENT << nl;
    #if (OPENFOAM >= 1812)
    realAmplitude_.writeEntry("realAmplitude", os);
    imagAmplitude_.writeEntry("imagAmplitude", os);
    this->writeEntry("value", os);
    #else
    writeEntry(os, "realAmplitude", realAmplitude_);
    writeEntry(os, "imagAmplitude", imagAmplitude_);
    writeEntry(os, "value", *this);
    #endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    complexAmplitudeDisplacementPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
