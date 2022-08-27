/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Jiri Furst, CTU in Prague
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
    Foam::functionObjects::mixingInterface

Group
    grpFieldFunctionObjects

Description
    Mixing interface implements a mixing plane boundary condition for MRF calculations.
    The condition connects two patches, namely an 'upstreamPatch' and 'downstreamPatch'
    and transfers chosen quantities (after mixing) in given direction.
    
    The condition is especially tailored for compressible flows and therefore it transfers:
    - total pressure:     downstream
    - total temperature:  downstream
    - velocity direction: downstream
    - static pressure:    upstream
    - other quantities:   downstream

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    mixingInterface1
    {
        // Mandatory entries (unmodifiable)
        type           mixingInterface;
        libs           (myFunctionObjects);


        upstreamPatch      statorOutlet;
        downstreamPatch    rotorInlet;
        axis               (1 0 0);     // Rotation axis of the machne
        configuration      axial;       // axial for axial machines, radial for radial machines
        segments           10;          // Splits patches into 10 segments
        convective         (k omega);   // Other convective scalar quantities
        relaxation         0.5;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                        | Type | Req'd | Dflt
      type       | Type name: mixingInterface         | word |  yes  | -
      libs       | Library name: mixingInterfaceFunctionObject  <!--
             --> | word  | yes | -
      upstreamPatch   | Name of the upstream patch    | word   |  yes  | -
      downstreamPatch | Name of the downstream patch  | word   |  yes  | -
      axis            | Rotation axis                 | vector |  yes  | -
      configuration   | Radial or axial turbine       | word   |  yes  | -
      segments        | Number of segments for mixing | label  |  yes  | 10
      convective      | Other convective quantities   | wordList |  no  | -
      relaxation      | Relaxation parameter          | scalar | no    | 1
      frequency       | Frequency of mixing update    | label  | no    | 1
    \endtable

    Options for the \c configuration entry:
    \verbatim
        axial   | for axial turbine, parametrization is done w.r. to radius
        radial  | for axial turbine, parametrization is done w.r. to axial coordinate
      ...
    \endverbatim

    The object assumes that the boundary conditions for set as:
    - upstreamPatch:
      - zeroGradient for U, T and other convective quantities
      - fixedValue for p
    - downstreamPatch:
      - totalPressure and totalTemperature for p and T
      - subsonicInletTotal for U (for myLusgsFoam)
      - fixedValue for other convective quantities
    Usage by the \c postProcess utility is not available.

Note
    - The object is fully parallelized

See also
    - Foam::functionObject

SourceFiles
    mixingInterface.C

\*---------------------------------------------------------------------------*/

#ifndef mixingInterface_H
#define mixingInterface_H

#include "regionFunctionObject.H"
#include "volFieldsFwd.H"
#include "vectorField.H"

#include <functional>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class mixingInterface Declaration
\*---------------------------------------------------------------------------*/

class mixingInterface
:
    public regionFunctionObject
{
    // Private Data

    //- Options for configuration
    enum configurationType
    {
        AXIAL, //!< "Quantities are expressed as functions of distance from axis"
        RADIAL   //!< "Quantities are expressed as functions of axial position"
    };

    //- Name of the upstreamPatch
    const word upstreamPatch_;

    //- Name of the downstreamPatch
    const word downstreamPatch_;

    //- Axis of rotation
    vector axis_;

    //- Origin of rotation
    vector origin_;

    //- Configuration (axial/radial)
    enum configurationType configuration_;

    //- Number of segments for mixing
    const label segments_;

    //- List of additional convective variables
    const wordList convectiveVariables_;

    //- Parametrization (returns radius for axial config and axial coord for radial config)
    std::function<scalar(vector)> parameter_;

    //- Relaxation parameter (new BC value = (1-relax)*OLD + relax*NEW)
    scalar relax_;

    //- Frequency of updates
    label frequency_;
    
    //- Parameter bins
    scalarList paramBins_;

    //- ID of the upstream patch
    const label upstreamPatchID_;
    
    //- ID of the upstream patch
    const label downstreamPatchID_;

    //- Iteration number
    label counter_;
    
    // Private member functions

    //- Cast reference of objectRegistry to fvMesh
    const fvMesh& mesh() const;

    //- Get averages for each bin
    List<scalar> getScalarAverages(
        const scalarField& values,
        const scalarField& weights,
        const vectorField& coords
    ) const;

    //- Get averages for each bin, vectors are transformed to (x, r, theta)
    List<vector> getVectorAverages(
        const vectorField& values,
        const scalarField& weights,
        const vectorField& coords
    ) const;

    scalar interpolateScalar(const List<scalar>& y, scalar p) const;

    vector interpolateVector(const List<vector>& y, scalar p) const;

    label getBin(scalar p) const;
    
public:

    //- Runtime type information
    TypeName("mixingInterface");


    // Constructors

        //- Construct from Time and dictionary
        mixingInterface
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        mixingInterface(const mixingInterface&) = delete;

        //- No copy assignment
        void operator=(const mixingInterface&) = delete;


    //- Destructor
    virtual ~mixingInterface() = default;


    // Member Functions

        //- Read the mixingInterface data
        virtual bool read(const dictionary& dict);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        virtual bool end();

        //- Write the mixingInterface
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
