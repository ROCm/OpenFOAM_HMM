/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2019 OpenCFD Ltd.
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
    std::reference_wrapper<const labelList> ids_;


public:

    // Constructors

        //- Construct null
        meshedSurfRef()
        :
            points_(std::cref<pointField>(pointField::null())),
            faces_(std::cref<faceList>(faceList::null())),
            ids_(std::cref<labelList>(labelList::null()))
        {}


        //- Construct from components
        meshedSurfRef
        (
            const pointField& pts,
            const faceList& fcs,
            const labelList& ids = labelList::null()
        )
        :
            points_(std::cref<pointField>(pts)),
            faces_(std::cref<faceList>(fcs)),
            ids_(std::cref<labelList>(ids))
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
            return ids_.get();
        }

        //- Remove all references by redirecting to null objects
        void clear()
        {
            points_ = std::cref<pointField>(pointField::null());
            faces_ = std::cref<faceList>(faceList::null());
            ids_ = std::cref<labelList>(labelList::null());
        }

        //- Reset components
        void reset
        (
            const pointField& pts,
            const faceList& fcs,
            const labelList& ids = labelList::null()
        )
        {
            points_ = std::cref<pointField>(pts);
            faces_ = std::cref<faceList>(fcs);
            ids_ = std::cref<labelList>(ids);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
