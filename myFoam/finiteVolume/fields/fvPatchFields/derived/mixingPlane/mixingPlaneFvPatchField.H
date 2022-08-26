/*---------------------------------------------------------------------------*\
    =========                 |
    \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
     \\    /   O peration     |
      \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

    Class
    Foam::mixingPlaneFvPatchField

    Group
    grpBoundaryConditions

    Description
    This boundary condition provides a simple mixing plane boundary condition
    for turbomachinery. The distribution of a quantity is approximated by a polynomial
    approximation using least square method.  

    \heading Patch usage

    \table
    Property   | Description                    | Required | Default value
    source     | name of the neighbour patch    | yes      |
    order      | order of the polynomial        | yes      | 
    origin     | point at the axis of rotation  | no       | (0 0 0)
    axis       | the axis of rotation           | no       | (1 0 0)   
    parametrization | radial or axial           | no       | radial 
    \endtable

    Example of the boundary condition specification for patch at rotor side of a
    stator-rotor interface
    \verbatim
    myPatch
    {
    type            mixingPlane;
    source          interface_stator_rotor;
    order           1;   // linear polynomial
    origin          (0 0 0);
    axis            (1 0 0);
    value           uniform 0;
    parametrization radial;
    }
    \endverbatim

    Notes
    - the stator side should use zeroGradient
    - radial parametrization assumes values as f(r) (usefull for axial flow)
    - axial parametrization assumes values as f(x) (usefull for radial flow)
 
    SeeAlso
    Foam::fixedValueFvPatchField

    SourceFiles
    mixingPlaneFvPatchField.C

    \*---------------------------------------------------------------------------*/

#ifndef mixingPlaneFvPatchField_H
#define mixingPlaneFvPatchField_H

#include "fixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class mixingPlaneFvPatchField Declaration
    \*---------------------------------------------------------------------------*/

template<class Type>
class mixingPlaneFvPatchField
:
    public fixedValueFvPatchField<Type>
{
private:
    // Private meber function
    Type toXRTheta(const Type& data, const vector& rDir) const;
    
    Type fromXRTheta(const Type& data, const vector& rDir) const;
    

protected:

    // Protected data

    //- Options for parametrization
    enum parametrizationType
    {
        RADIAL, //!< "Quantities are expressed as functions of distance from axis"
        AXIAL   //!< "Quantities are expressed as functions of axial position"
    };

    //- origin (a point at the axis of rotation)
    const vector origin_;
    
    //- axis of rotation
    vector axis_;

    //- order of polynomial approximation
    const label order_;

    //- Name of the neighbour patch
    const word source_;

    //- Parametrization (axial/radial)
    enum parametrizationType parametrization_;

public:

    //- Runtime type information
    TypeName("mixingPlane");
    
    
    // Constructors
        
    //- Construct from patch and internal field 
    mixingPlaneFvPatchField
    (
        const fvPatch&,
        const DimensionedField<Type, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    mixingPlaneFvPatchField
    (
        const fvPatch&,
        const DimensionedField<Type, volMesh>&,
        const dictionary&
    );
    
    //- Construct by mapping given mixingPlaneFvPatchField onto a new patch
    mixingPlaneFvPatchField
    (
        const mixingPlaneFvPatchField<Type>&,
        const fvPatch&,
        const DimensionedField<Type, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Construct as copy
    mixingPlaneFvPatchField
    (
        const mixingPlaneFvPatchField<Type>&
    );

    //- Construct and return a clone
    virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
                (
                    new mixingPlaneFvPatchField<Type>(*this)
                );
        }

    //- Construct as copy setting internal field reference
    mixingPlaneFvPatchField
    (
        const mixingPlaneFvPatchField<Type>&,
        const DimensionedField<Type, volMesh>&
    );

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchField<Type> > clone
    (
        const DimensionedField<Type, volMesh>& iF
    ) const
        {
            return tmp<fvPatchField<Type> >
                (
                    new mixingPlaneFvPatchField<Type>(*this, iF)
                );
        }


    // Member functions

    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();

    //- Write
    virtual void write(Ostream&) const;

};

template<>
scalar mixingPlaneFvPatchField<scalar>::toXRTheta(const scalar& data, const vector& rDir) const;

template<>
scalar mixingPlaneFvPatchField<scalar>::fromXRTheta(const scalar& data, const vector& rDir) const;

template<>
vector mixingPlaneFvPatchField<vector>::toXRTheta(const vector& data, const vector& rDir) const;

template<>
vector mixingPlaneFvPatchField<vector>::fromXRTheta(const vector& data, const vector& rDir) const;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "mixingPlaneFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
