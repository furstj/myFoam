/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "gasProperties.H"
#include "HashTable.H"
#include "Switch.H"


#include <typeinfo>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gasProperties, 0);
    defineRunTimeSelectionTable(gasProperties, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


namespace Foam
{


gasProperties::~gasProperties()
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //


autoPtr<gasProperties> gasProperties::New
(
    const dictionary& dict
)
{
    if (debug)
    {
        InfoInFunction << "Constructing gasProperties" << endl;
    }


    const dictionary& thermoDict = dict.subDict("thermoType");
    
    const word transportName(thermoDict.lookup("transport"));
    const word thermoName(thermoDict.lookup("thermo"));
    const word EOSName(thermoDict.lookup("equationOfState"));
    const word specieName(thermoDict.lookup("specie"));
    const word energyName(thermoDict.lookup("energy"));

    const word gasPropertiesTypeName = "gasProperties<" + transportName + "<"
        + thermoName + "<" + EOSName + "<" + specieName + ">>," + energyName + ">>";
    
    auto cstrIter = dictionaryConstructorTablePtr_->cfind(gasPropertiesTypeName);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown gasProperties type "
                << gasPropertiesTypeName << nl << nl
                << "Valid gasProperties types are:" << nl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalError);
    }
    
    return autoPtr<gasProperties>(cstrIter()(dict));
}
    

}


// ************************************************************************* //
