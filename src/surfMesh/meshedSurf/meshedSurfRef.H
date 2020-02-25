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
    Foam::meshedSurfRef

Description
    Implements a meshed surface by referencing existing faces and points.

\*---------------------------------------------------------------------------*/

#ifndef meshedSurfRef_H
#define meshedSurfRef_H

#include "meshedSurf.H"
#include <functional>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class meshedSurfRef Declaration
\*---------------------------------------------------------------------------*/

class meshedSurfRef
:
    public meshedSurf
{
    std::reference_wrapper<const pointField> points_;
    std::reference_wrapper<const faceList> faces_;
    std::reference_wrapper<const labelList> zoneIds_;
    std::reference_wrapper<const labelList> faceIds_;


public:

    // Constructors

        //- Default construct
        meshedSurfRef()
        :
            points_(std::cref<pointField>(pointField::null())),
            faces_(std::cref<faceList>(faceList::null())),
            zoneIds_(std::cref<labelList>(labelList::null())),
            faceIds_(std::cref<labelList>(labelList::null()))
        {}


        //- Construct from components
        meshedSurfRef
        (
            const pointField& pointLst,
            const faceList& faceLst,
            const labelList& zoneIdLst = labelList::null(),
            const labelList& faceIdLst = labelList::null()
        )
        :
            points_(std::cref<pointField>(pointLst)),
            faces_(std::cref<faceList>(faceLst)),
            zoneIds_(std::cref<labelList>(zoneIdLst)),
            faceIds_(std::cref<labelList>(faceIdLst))
        {}


    //- Destructor
    virtual ~meshedSurfRef() = default;


    // Member Functions

        //- The points used for the surface
        virtual const pointField& points() const
        {
            return points_.get();
        }

        //- The faces used for the surface
        virtual const faceList& faces() const
        {
            return faces_.get();
        }

        //- Per-face zone/region information.
        virtual const labelList& zoneIds() const
        {
            return zoneIds_.get();
        }

        //- Per-face identifier (eg, element Id)
        virtual const labelList& faceIds() const
        {
            return faceIds_.get();
        }

        //- Remove all references by redirecting to null objects
        void clear()
        {
            points_ = std::cref<pointField>(pointField::null());
            faces_ = std::cref<faceList>(faceList::null());
            zoneIds_ = std::cref<labelList>(labelList::null());
            faceIds_ = std::cref<labelList>(labelList::null());
        }

        //- Reset components
        void reset
        (
            const pointField& pointLst,
            const faceList& faceLst,
            const labelList& zoneIdLst = labelList::null(),
            const labelList& faceIdLst = labelList::null()
        )
        {
            points_ = std::cref<pointField>(pointLst);
            faces_ = std::cref<faceList>(faceLst);
            zoneIds_ = std::cref<labelList>(zoneIdLst);
            faceIds_ = std::cref<labelList>(faceIdLst);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
