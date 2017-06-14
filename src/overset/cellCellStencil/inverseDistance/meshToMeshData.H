/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::meshToMeshData

Description
    Transfers refinement levels such that slow transition between levels is
    maintained. Used in FaceCellWave.

SourceFiles
    meshToMeshDataI.H
    meshToMeshData.C

\*---------------------------------------------------------------------------*/

#ifndef meshToMeshData_H
#define meshToMeshData_H

#include "point.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyPatch;
class polyMesh;


// Forward declaration of friend functions and operators

class meshToMeshData;

Istream& operator>>(Istream&, meshToMeshData&);
Ostream& operator<<(Ostream&, const meshToMeshData&);


/*---------------------------------------------------------------------------*\
                           Class meshToMeshData Declaration
\*---------------------------------------------------------------------------*/

class meshToMeshData
{
public:

    //- Class used to pass non-cell data to the update function
    class trackData
    {
    public:

        const polyMesh& tgtMesh_;


        // Constructors

            trackData(const polyMesh& tgtMesh)
            :
                tgtMesh_(tgtMesh)
            {}

    };


    // Private data

        // Corresponding cell in tgtMesh
        label tgtCelli_;


public:

    // Constructors

        //- Construct null
        inline meshToMeshData();

        //- Construct from count
        inline meshToMeshData
        (
            const label tgtCelli
        );


    // Member Functions

        // Access

            inline label tgtCell() const
            {
                return tgtCelli_;
            }

            inline label& tgtCell()
            {
                return tgtCelli_;
            }


        // Needed by FaceCellWave

            //- Check whether origin has been changed at all or
            //  still contains original (invalid) value.
            template<class TrackingData>
            inline bool valid(TrackingData&) const;

            //- Check for identical geometrical data. Used for cyclics checking.
            template<class TrackingData>
            inline bool sameGeometry
            (
                const polyMesh&,
                const meshToMeshData&,
                const scalar,
                TrackingData&
            ) const;

            //- Convert any absolute coordinates into relative to (patch)face
            //  centre
            template<class TrackingData>
            inline void leaveDomain
            (
                const polyMesh&,
                const polyPatch&,
                const label patchFacei,
                const point& faceCentre,
                TrackingData&
            );

            //- Reverse of leaveDomain
            template<class TrackingData>
            inline void enterDomain
            (
                const polyMesh&,
                const polyPatch&,
                const label patchFacei,
                const point& faceCentre,
                TrackingData&
            );

            //- Apply rotation matrix to any coordinates
            template<class TrackingData>
            inline void transform
            (
                const polyMesh&,
                const tensor&,
                TrackingData&
            );

            //- Influence of neighbouring face.
            template<class TrackingData>
            inline bool updateCell
            (
                const polyMesh&,
                const label thisCelli,
                const label neighbourFacei,
                const meshToMeshData& neighbourInfo,
                const scalar tol,
                TrackingData&
            );

            //- Influence of neighbouring cell.
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const label neighbourCelli,
                const meshToMeshData& neighbourInfo,
                const scalar tol,
                TrackingData&
            );

            //- Influence of different value on same face.
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const meshToMeshData& neighbourInfo,
                const scalar tol,
                TrackingData&
            );

            //- Same (like operator==)
            template<class TrackingData>
            inline bool equal
            (
                const meshToMeshData&,
                TrackingData&
            ) const;

    // Member Operators

        // Needed for List IO
        inline bool operator==(const meshToMeshData&) const;

        inline bool operator!=(const meshToMeshData&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const meshToMeshData&);
        friend Istream& operator>>(Istream&, meshToMeshData&);
};


//- Data associated with meshToMeshData type are contiguous
template<>
inline bool contiguous<meshToMeshData>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshToMeshDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
