/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    Foam::ensightPartNonMeshFaces

Description
    An implementation of ensightPart to work on self-contained faces and points
    (without a mesh).

SourceFiles
    ensightPartNonMeshFaces.C

\*---------------------------------------------------------------------------*/

#ifndef ensightPartNonMeshFaces_H
#define ensightPartNonMeshFaces_H

#include "ensightPartFaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ensightPartNonMeshFaces Declaration
\*---------------------------------------------------------------------------*/

class ensightPartNonMeshFaces
:
    public ensightPartFaces
{
    // Private data

        //- faces (reference)
        const faceList& faces_;

        //- points (reference)
        const pointField& points_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const ensightPartNonMeshFaces&);

        //- track points used
        virtual localPoints calcLocalPoints() const;

public:

    //- Runtime type information
    TypeName("ensightNonMeshFaces");

    // Constructors

        //- Construct from faces and points
        ensightPartNonMeshFaces
        (
            label partNumber,
            const string& partDescription,
            const faceList& faces,
            const pointField& points
        );

        //- Construct as copy
        ensightPartNonMeshFaces(const ensightPartNonMeshFaces& part);

        //- Construct from Istream
        ensightPartNonMeshFaces(Istream& is);

        //- Construct on freestore from Istream
        static autoPtr<ensightPartNonMeshFaces> New(Istream& is)
        {
            return autoPtr<ensightPartNonMeshFaces>
            (
                new ensightPartNonMeshFaces(is)
            );
        }


    //- Destructor
    virtual ~ensightPartNonMeshFaces();


    // Member Functions

        //- element connectivity
        virtual void writeConnectivity
        (
            ensightGeoFile& os,
            const word& key,
            const labelList& idList,
            const labelList& pointMap
        ) const;

        virtual void writeGeometry(ensightGeoFile& os) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
