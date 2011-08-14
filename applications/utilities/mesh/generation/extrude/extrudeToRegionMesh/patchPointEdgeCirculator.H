/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::patchPointEdgeCirculator

Description
    Walks from starting edge/face around point on patch.

    -# Use in-place: \n
        \code
            patchPointEdgeCirculator circ(..);

            // Walk
            do
            {
                Info<< "edge:" << circ.edgeID() << endl;
                ++circ;
            }
            while (circ != circ.end());
        \endcode

    -# Use like STL iterator: \n
        \code
            patchPointEdgeCirculator circ(..);
            for
            (
                patchPointEdgeCirculator iter = circ.begin();
                iter != circ.end();
                ++iter
            )
            {
                Info<< "edge:" << iter.edgeID() << endl;
            }
        \endcode


SourceFiles
    patchPointEdgeCirculator.C

\*---------------------------------------------------------------------------*/

#ifndef patchPointEdgeCirculator_H
#define patchPointEdgeCirculator_H

#include "face.H"
#include "ListOps.H"
#include "primitiveFacePatch.H"
#include "PackedBoolList.H"
#include "InfoProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class patchPointEdgeCirculator Declaration
\*---------------------------------------------------------------------------*/

class patchPointEdgeCirculator
{
    // Static data members

        //- end iterator
        static const patchPointEdgeCirculator endConstIter;


    // Private data

        //- Patch
        const primitiveFacePatch& patch_;

        const PackedBoolList& nonManifoldEdge_;

        //- Current edge
        label edgeID_;

        //- Current direction (face, expressed as index into edgeFaces()[edgeID]
        label index_;

        //- Point
        label pointID_;

        //- Starting edge so we know when to stop.
        label startEdgeID_;


    // Private Member Functions

        //- Set to end() iterator
        inline void setEnd();

        //- Set edgeID_ to be the other edge on the face that uses the
        //  point.
        inline void otherEdge(const label cellI);

public:

    // Constructors

        //- Construct from components
        inline patchPointEdgeCirculator
        (
            const primitiveFacePatch&,
            const PackedBoolList& nonManifoldEdge,
            const label edgeID,
            const label index,
            const label pointID
        );

        //- Construct as copy
        inline patchPointEdgeCirculator(const patchPointEdgeCirculator&);


    // Member Functions

        inline label edgeID() const;

        inline label index() const;

        inline label pointID() const;

        //- Helper: get face according to the index.
        //  Returns -1 if on end.
        inline label faceID() const;

        //- Walk back until non-manifold edge (if any) or minimum edge.
        inline void setCanonical();


    // Member Operators

        inline void operator=(const patchPointEdgeCirculator& iter);

        inline bool operator==(const patchPointEdgeCirculator& iter) const;

        inline bool operator!=(const patchPointEdgeCirculator& iter) const;

        //- Step to next face.
        inline patchPointEdgeCirculator& operator++();

        //- iterator set to the beginning face. For internal edges this is
        //  the current face. For boundary edges this is the first boundary face
        //  reached from walking back (i.e. in opposite direction to ++)
        inline patchPointEdgeCirculator begin() const;
        inline patchPointEdgeCirculator cbegin() const;

        //- iterator set to beyond the end of the walk.
        inline const patchPointEdgeCirculator& end() const;
        inline const patchPointEdgeCirculator& cend() const;

    // Info

        //- Return info proxy.
        //  Used to print token information to a stream
        InfoProxy<patchPointEdgeCirculator> info() const
        {
            return *this;
        }

        friend Ostream& operator<<
        (
            Ostream&,
            const InfoProxy<patchPointEdgeCirculator>&
        );
};

Ostream& operator<<(Ostream&, const InfoProxy<patchPointEdgeCirculator>&);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "patchPointEdgeCirculatorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
