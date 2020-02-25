/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::mergedSurf

Description
    Simple class to manage surface merging information.

    Merging is done with PatchTools::gatherAndMerge()

SourceFiles
    mergedSurf.C

\*---------------------------------------------------------------------------*/

#ifndef mergedSurf_H
#define mergedSurf_H

#include "meshedSurf.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class mergedSurf Declaration
\*---------------------------------------------------------------------------*/

class mergedSurf
:
    public meshedSurf
{
    pointField points_;
    faceList   faces_;
    labelList  pointsMap_;

    labelList  zoneIds_;
    labelList  faceIds_;


public:

    // Constructors

        //- Default construct
        mergedSurf() = default;

        //- Copy construct
        mergedSurf(const mergedSurf&) = default;

        //- Move construct
        mergedSurf(mergedSurf&&) = default;

        //- Construct and merge
        mergedSurf
        (
            const meshedSurf& unmergedSurface,
            const scalar mergeDim
        );

        //- Construct and merge
        mergedSurf
        (
            const pointField& unmergedPoints,
            const faceList& unmergedFaces,
            const scalar mergeDim
        );

        //- Construct and merge
        mergedSurf
        (
            const pointField& unmergedPoints,
            const faceList& unmergedFaces,
            const labelList& origZoneIds,
            const labelList& origFaceIds,
            const scalar mergeDim
        );


    //- Destructor
    virtual ~mergedSurf() = default;


    // Access Member Functions

        //- Can use (parallel only)
        static bool use();

        //- Number of faces
        label size() const
        {
            return faces_.size();
        }

        //- Const access to (global) points used for the surface
        virtual const pointField& points() const
        {
            return points_;
        }

        //- Const access to the surface faces
        virtual const faceList& faces() const
        {
            return faces_;
        }

        //- Per-face zone/region information
        virtual const labelList& zoneIds() const
        {
            return zoneIds_;
        }

        //- Per-face identifier (eg, element Id)
        virtual const labelList& faceIds() const
        {
            return faceIds_;
        }

        //- Map for reordered points (old-to-new)
        const labelList& pointsMap() const
        {
            return pointsMap_;
        }


    // Edit

        //- Clear all storage
        void clear();

        //- Merge meshed surfaces (in parallel only).
        bool merge
        (
            const meshedSurf& unmergedSurface,
            const scalar mergeDim
        );

        //- Merge meshed surfaces (in parallel only).
        bool merge
        (
            const pointField& unmergedPoints,
            const faceList& unmergedFaces,
            const scalar mergeDim
        );

        //- Merge meshed surfaces (in parallel only).
        bool merge
        (
            const pointField& unmergedPoints,
            const faceList& unmergedFaces,
            const labelList& origZoneIds,
            const labelList& origFaceIds,
            const scalar mergeDim
        );


    // Member Operators

        //- Copy assignment
        mergedSurf& operator=(const mergedSurf&) = default;

        //- Move assignment
        mergedSurf& operator=(mergedSurf&&) = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
