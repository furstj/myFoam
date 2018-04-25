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
    Foam::PointData

Description
    Variant of pointEdgePoint with some transported additional data. Templated
    on the transported data type.

SourceFiles
    PointDataI.H
    PointData.C

\*---------------------------------------------------------------------------*/

#ifndef PointData_H
#define PointData_H

#include "pointEdgePoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;
template<class DataType>
class PointData;

// Forward declaration of friend functions and operators
template<class DataType>
Ostream& operator<<(Ostream&, const PointData<DataType>&);
template<class DataType>
Istream& operator>>(Istream&, PointData<DataType>&);

/*---------------------------------------------------------------------------*\
                           Class PointData Declaration
\*---------------------------------------------------------------------------*/

template<class DataType>
class PointData
:
    public pointEdgePoint
{
private:

    // Private data

        //- Additional transported data
        DataType data_;


public:

    // Constructors

        //- Construct null
        inline PointData();

        //- Construct from origin, distance and data
        inline PointData
        (
            const point& origin,
            const scalar distSqr,
            const DataType& data
        );

        //- Construct as copy
        inline PointData(const PointData&);


    // Member Functions

        // Access

            //- Const access the data
            inline const DataType& data() const;


        // Needed by meshWave

            //- Apply rotation matrix to the data
            template<class TrackingData>
            inline void transform
            (
                const tensor& rotTensor,
                TrackingData& td
            );

            //- Influence of edge on point
            template<class TrackingData>
            inline bool updatePoint
            (
                const polyMesh& mesh,
                const label pointI,
                const label edgeI,
                const PointData<DataType>& edgeInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same point.
            //  Merge new and old info.
            template<class TrackingData>
            inline bool updatePoint
            (
                const polyMesh& mesh,
                const label pointI,
                const PointData<DataType>& newPointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same point.
            //  No information about current position whatsoever.
            template<class TrackingData>
            inline bool updatePoint
            (
                const PointData<DataType>& newPointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of point on edge.
            template<class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const label edgeI,
                const label pointI,
                const PointData<DataType>& pointInfo,
                const scalar tol,
                TrackingData& td
            );

    // Member Operators

        inline bool operator==(const PointData<DataType>&) const;
        inline bool operator!=(const PointData<DataType>&) const;


    // IOstream Operators

        friend Ostream& operator<< <DataType>
        (
            Ostream&,
            const PointData<DataType>&
        );
        friend Istream& operator>> <DataType>
        (
            Istream&,
            PointData<DataType>&
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Data associated with PointData types is contiguous

template<>
inline bool contiguous<PointData<scalar>>()
{
    return true;
}

template<>
inline bool contiguous<PointData<vector>>()
{
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PointDataI.H"

#ifdef NoRepository
#   include "PointData.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
