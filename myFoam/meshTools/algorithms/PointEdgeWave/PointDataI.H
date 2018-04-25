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

#include "polyMesh.H"
#include "transform.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class DataType>
inline Foam::PointData<DataType>::PointData()
:
    pointEdgePoint()
{}


template<class DataType>
inline Foam::PointData<DataType>::PointData
(
    const point& origin,
    const scalar distSqr,
    const DataType& data
)
:
    pointEdgePoint(origin, distSqr),
    data_(data)
{}


template<class DataType>
inline Foam::PointData<DataType>::PointData(const PointData<DataType>& wpt)
:
    pointEdgePoint(wpt),
    data_(wpt.data())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class DataType>
inline const DataType& Foam::PointData<DataType>::data() const
{
    return data_;
}


template<class DataType>
template<class TrackingData>
inline void Foam::PointData<DataType>::transform
(
    const tensor& rotTensor,
    TrackingData& td
)
{
    pointEdgePoint::transform(rotTensor, td);
    data_ = Foam::transform(rotTensor, data_);
}


template<class DataType>
template<class TrackingData>
inline bool Foam::PointData<DataType>::updatePoint
(
    const polyMesh& mesh,
    const label pointI,
    const label edgeI,
    const PointData<DataType>& edgeInfo,
    const scalar tol,
    TrackingData& td
)
{
    if
    (
        pointEdgePoint::updatePoint
        (
            mesh,
            pointI,
            edgeI,
            edgeInfo,
            tol,
            td
        )
    )
    {
        data_ = edgeInfo.data_;

        return true;
    }
    else
    {
        return false;
    }
}


template<class DataType>
template<class TrackingData>
inline bool Foam::PointData<DataType>::updatePoint
(
    const polyMesh& mesh,
    const label pointI,
    const PointData<DataType>& newPointInfo,
    const scalar tol,
    TrackingData& td
)
{
    if
    (
        pointEdgePoint::updatePoint
        (
            mesh,
            pointI,
            newPointInfo,
            tol,
            td
        )
    )
    {
        data_ = newPointInfo.data_;

        return true;
    }
    else
    {
        return false;
    }
}


template<class DataType>
template<class TrackingData>
inline bool Foam::PointData<DataType>::updatePoint
(
    const PointData<DataType>& newPointInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (pointEdgePoint::updatePoint(newPointInfo, tol, td))
    {
        data_ = newPointInfo.data_;

        return true;
    }
    else
    {
        return false;
    }
}


template<class DataType>
template<class TrackingData>
inline bool Foam::PointData<DataType>::updateEdge
(
    const polyMesh& mesh,
    const label edgeI,
    const label pointI,
    const PointData<DataType>& pointInfo,
    const scalar tol,
    TrackingData& td

)
{
    if
    (
        pointEdgePoint::updateEdge
        (
            mesh,
            edgeI,
            pointI,
            pointInfo,
            tol,
            td
        )
    )
    {
        data_ = pointInfo.data_;

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class DataType>
inline bool Foam::PointData<DataType>::operator==
(
    const Foam::PointData<DataType>& rhs
)
const
{
    return pointEdgePoint::operator==(rhs) && (data() == rhs.data());
}


template<class DataType>
inline bool Foam::PointData<DataType>::operator!=
(
    const Foam::PointData<DataType>& rhs
)
const
{
    return !(*this == rhs);
}


// ************************************************************************* //
