/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::PointIntegrateData

Description
    Integrate along selected edges using PointEdgeWave.

SourceFiles
    PointIntegrateDataI.H

\*---------------------------------------------------------------------------*/

#ifndef PointIntegrateData_H
#define PointIntegrateData_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;
template<class DataType>
class PointIntegrateData;

// Forward declaration of friend functions and operators
template<class DataType>
Ostream& operator<<(Ostream&, const PointIntegrateData<DataType>&);
template<class DataType>
Istream& operator>>(Istream&, PointIntegrateData<DataType>&);

/*---------------------------------------------------------------------------*\
                     Class PointIntegrateData declaration
\*---------------------------------------------------------------------------*/

template<class DataType>
class PointIntegrateData
{
private:

    // Private data

        //- Valid flag
        bool valid_;

        //- Integrated data
        DataType data_;


public:

    //- Class used to pass extra data
    class trackingData
    {
    public:

        UList<DataType>& edgeData_;

        trackingData(UList<DataType>& edgeData)
        :
            edgeData_(edgeData)
        {}
    };



    // Constructors

        //- Construct null
        inline PointIntegrateData();

        //- Construct from data
        inline PointIntegrateData(const DataType& data);


    // Member Functions

        // Access

            //- Const access the data
            inline const DataType& data() const
            {
                return data_;
            };


        // Needed by PointEdgeWave

            //- Check whether original (invalid) value.
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Check for identical geometrical data. Used for cyclics checking.
            template<class TrackingData>
            inline bool sameGeometry
            (
                const PointIntegrateData<DataType>&,
                const scalar tol,
                TrackingData& td
            ) const;

            //- Convert origin to relative vector to leaving point
            //  (= point coordinate)
            template<class TrackingData>
            inline void leaveDomain
            (
                const polyPatch& patch,
                const label patchPointi,
                const point& pos,
                TrackingData& td
            );

            //- Convert relative origin to absolute by adding entering point
            template<class TrackingData>
            inline void enterDomain
            (
                const polyPatch& patch,
                const label patchPointi,
                const point& pos,
                TrackingData& td
            );

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
                const PointIntegrateData<DataType>& edgeInfo,
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
                const PointIntegrateData<DataType>& newPointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same point.
            //  No information about current position whatsoever.
            template<class TrackingData>
            inline bool updatePoint
            (
                const PointIntegrateData<DataType>& newPointInfo,
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
                const PointIntegrateData<DataType>& pointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Same (like operator==)
            template<class TrackingData>
            inline bool equal
            (
                const PointIntegrateData<DataType>&,
                TrackingData& td
            ) const;


    // Member Operators

        // Needed for List IO
        inline bool operator==(const PointIntegrateData<DataType>&) const;
        inline bool operator!=(const PointIntegrateData<DataType>&) const;


    // IOstream Operators

        friend Ostream& operator<< <DataType>
        (
            Ostream&,
            const PointIntegrateData<DataType>&
        );
        friend Istream& operator>> <DataType>
        (
            Istream&,
            PointIntegrateData<DataType>&
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Data associated with PointIntegrateData types is contiguous

template<>
inline bool contiguous<PointIntegrateData<scalar>>()
{
    return true;
}

template<>
inline bool contiguous<PointIntegrateData<vector>>()
{
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PointIntegrateDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
