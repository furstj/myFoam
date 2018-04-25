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

\*---------------------------------------------------------------------------*/

#include "PointData.H"

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class DataType>
Foam::Ostream& Foam::operator<<(Ostream& os, const PointData<DataType>& pd)
{
    if (os.format() == IOstream::ASCII)
    {
        return os
            << static_cast<const pointEdgePoint&>(pd)
            << token::SPACE << pd.data();
    }
    else
    {
        return os
            << static_cast<const pointEdgePoint&>(pd)
            << pd.data();
    }
}


template<class DataType>
Foam::Istream& Foam::operator>>(Istream& is, PointData<DataType>& pd)
{
    return is >> static_cast<pointEdgePoint&>(pd) >> pd.data_;
}


// ************************************************************************* //
