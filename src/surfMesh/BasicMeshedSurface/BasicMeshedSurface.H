/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    Foam::BasicMeshedSurface

Description
    Holds surfaces without any zone information

SourceFiles
    BasicMeshedSurface.C

\*---------------------------------------------------------------------------*/

#ifndef BasicMeshedSurface_H
#define BasicMeshedSurface_H

#include "PrimitivePatch.H"
#include "PrimitivePatchTools.H"
#include "pointField.H"
#include "face.H"
#include "triFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

/*---------------------------------------------------------------------------*\
                     Class BasicMeshedSurface Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class BasicMeshedSurface
:
    public PrimitivePatch<Face, ::Foam::List, pointField, point>
{

    //- Typedefs for convenience
        typedef PrimitivePatch
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

        //- Set new regions from faceMap
        virtual void remapFaces(const UList<label>& faceMap);

public:

    // Static

        //- Face storage only handles triangulated faces
        inline static bool isTri();


    // Constructors

        //- Construct null
        BasicMeshedSurface();

        //- Construct by transferring components (points, faces).
        BasicMeshedSurface
        (
            const Xfer<pointField>&,
            const Xfer<List<Face> >&
        );

    // Destructor

        virtual ~BasicMeshedSurface();

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
            const Xfer<pointField>&,
            const Xfer<List<Face> >&
        );

        //- Transfer components (points, faces).
        virtual void reset
        (
            const Xfer<List<point> >&,
            const Xfer<List<Face> >&
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
        //  Optionally returning a map of original face Ids (zero-sized when
        //  no triangulation was done)
        virtual label triangulate();

        //- Triangulate in-place, setting a map of original face Ids
        //  faceMap is zero-sized when no triangulation was done
        virtual label triangulate(List<label>& faceMap);


    // Write

        void writeStats(Ostream& os) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Specialization for holding triangulated information
template<>
inline bool BasicMeshedSurface<triFace>::isTri()
{
    return true;
}

//- Specialization for holding triangulated information
template<>
inline label BasicMeshedSurface<triFace>::triangulate()
{
    return 0;
}

//- Specialization for holding triangulated information
template<>
inline label BasicMeshedSurface<triFace>::triangulate(List<label>& faceMap)
{
    if (&faceMap)
    {
        faceMap.clear();
    }

    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BasicMeshedSurface.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
