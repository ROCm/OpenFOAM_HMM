/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::mergedSurf

Description
    Simple class to manage surface merging information

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
    labelList  zones_;
    labelList  pointsMap_;

    //- Disallow default bitwise copy construct
    mergedSurf(const mergedSurf&) = delete;

    // Assignment is needed for lists

public:

    // Constructors

        //- Construct null
        mergedSurf();

        //- Construct and merge meshed surfaces immediately (in parallel only).
        mergedSurf(const meshedSurf&, const scalar mergeDim);


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

        //- Const access to per-face zone/region information
        virtual const labelList& zoneIds() const
        {
            return zones_;
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
        bool merge(const meshedSurf&, const scalar mergeDim);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
