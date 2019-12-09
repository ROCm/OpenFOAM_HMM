/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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
    Foam::ensightPartFaces

Description
    An implementation of ensightPart to hold mesh faces.

SourceFiles
    ensightPartFaces.C

\*---------------------------------------------------------------------------*/

#ifndef ensightPartFaces_H
#define ensightPartFaces_H

#include "ensightPart.H"
#include "ensightFaces.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ensightPartFaces Declaration
\*---------------------------------------------------------------------------*/

class ensightPartFaces
:
    public ensightFaces,
    public ensightPart
{
    // Private data

        //- Start offset for patch
        const label start_;

        //- Patch index
        const label patchIndex_;

        //- The referenced faces
        const faceList& faces_;

        //- The referenced pointField
        const pointField& points_;

        //- Can skip local point renumbering when points are contiguous
        const bool contiguousPoints_;


    // Private Member Functions

        //- Track points used
        localPoints calcLocalPoints() const;

        //- Element connectivity
        void writeConnectivity
        (
            ensightGeoFile&,
            const word& key,
            const labelUList& idList,
            const labelUList& pointMap
        ) const;


        //- Helper: write connectivity
        void writeConnectivity
        (
            ensightGeoFile&,
            const word& key,
            const faceList&,
            const labelUList& idList,
            const labelUList& pointMap
        ) const;


        //- No copy construct
        ensightPartFaces(const ensightPartFaces&) = delete;

        //- No copy assignment
        void operator=(const ensightPartFaces&) = delete;


public:

    //- Runtime type information
    TypeName("ensightFaces");


    // Constructors

        //- Construct part with 0-based index, description, points and faces
        //  Can skip local point renumbering when points are contiguous
        ensightPartFaces
        (
            label partIndex,
            const string& description,
            const pointField& points,
            const faceList& faces,
            const bool contiguousPoints = false
        );

        //- Construct from polyMesh and polyPatch
        //- Part receives the name of the patch unless otherwise specified.
        ensightPartFaces
        (
            label partIndex,
            const polyMesh& mesh,
            const polyPatch& patch,
            const string& partName = ""
        );

        //- Construct from polyPatch
        //- Part receives the name of the patch unless otherwise specified.
        ensightPartFaces
        (
            label partIndex,
            const polyPatch& patch,
            const string& partName = ""
        );


    //- Destructor
    virtual ~ensightPartFaces() = default;


    // Member Functions

    // Access

        //- Part index (0-based)
        virtual label index() const
        {
            return ensightFaces::index();
        }


        //- Number of elements in this part
        virtual label size() const
        {
            return ensightFaces::size();
        }


        //- Return the patch index, -1 when not in use.
        inline label patchIndex() const
        {
            return patchIndex_;
        }


    // Output

        //- Write summary information about the object
        virtual void writeSummary(Ostream& os) const;

        //- Print various types of debugging information
        virtual void dumpInfo(Ostream& os) const;

        //- Write geometry
        virtual void write(ensightGeoFile& os) const;

        //- Helper: write geometry given the pointField
        virtual void write(ensightGeoFile& os, const pointField& points) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
