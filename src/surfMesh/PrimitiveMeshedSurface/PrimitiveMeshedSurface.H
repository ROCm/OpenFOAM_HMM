/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::PrimitiveMeshedSurface

Description
    Holds surfaces without any patch information

SourceFiles
    PrimitiveMeshedSurface.C

\*---------------------------------------------------------------------------*/

#ifndef PrimitiveMeshedSurface_H
#define PrimitiveMeshedSurface_H

#include "PrimitivePatchExtra.H"
#include "pointField.H"
#include "xfer.H"
#include "face.H"
#include "triFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

/*---------------------------------------------------------------------------*\
                  Class PrimitiveMeshedSurface Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class PrimitiveMeshedSurface
:
    public PrimitivePatchExtra<Face, ::Foam::List, pointField, point>
{

    //- Typedefs for convenience
        typedef PrimitivePatchExtra
        <
            Face,
            ::Foam::List,
            pointField,
            point
        >
        ParentType;

protected:
    // Protected Member Functions

        //- Non-const access to global points
        pointField& storedPoints()
        {
            return const_cast<pointField&>(ParentType::points());
        }

        //- Non-const access to the faces
        List<Face>& storedFaces()
        {
            return static_cast<List<Face> &>(*this);
        }

        //- Set new regions/patches from faceMap
        virtual void remapRegions(List<label>& faceMap);

public:

    // Static

        //- Face storage only handles triangulated faces
        inline static bool isTri();

    // Constructors

        //- Construct null
        PrimitiveMeshedSurface();

        //- Construct by transferring components (points, faces).
        PrimitiveMeshedSurface
        (
            const xfer<pointField>&,
            const xfer<List<Face> >&
        );

    // Destructor

        virtual ~PrimitiveMeshedSurface();

    // Member Functions

    // Access

        //- Return const access to the faces
        inline const List<Face>& faces() const
        {
            return static_cast<const List<Face> &>(*this);
        }

    // Edit

        //- Clear all storage
        virtual void clear();

        //- Move points
        virtual void movePoints(const pointField&);

        //- Scale points. A non-positive factor is ignored
        virtual void scalePoints(const scalar&);

        //- Transfer components (points, faces).
        virtual void reset
        (
            const xfer<pointField>&,
            const xfer<List<Face> >&
        );

        //- Remove invalid faces
        virtual void cleanup(const bool verbose);

        virtual bool stitchFaces
        (
            const scalar tol=SMALL,
            const bool verbose=false
        );

        virtual bool checkFaces
        (
            const bool verbose=false
        );

        //- Triangulate in-place
        //  Returning the number of triangles added
        virtual label triangulate();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Specialization for holding triangulated information
template<>
inline bool PrimitiveMeshedSurface<triFace>::isTri()
{
    return true;
}

//- Specialization for holding triangulated information
template<>
inline label PrimitiveMeshedSurface<triFace>::triangulate()
{
    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PrimitiveMeshedSurface.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
