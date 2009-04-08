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
    Foam::featureEdgeMesh

Description

SourceFiles
    featureEdgeMeshI.H
    featureEdgeMesh.C

\*---------------------------------------------------------------------------*/

#ifndef featureEdgeMesh_H
#define featureEdgeMesh_H

#include "edgeMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class featureEdgeMesh Declaration
\*---------------------------------------------------------------------------*/

class featureEdgeMesh
:
    public edgeMesh
{
    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        featureEdgeMesh(const featureEdgeMesh&);

        //- Disallow default bitwise assignment
        void operator=(const featureEdgeMesh&);


public:

    // Constructors

        //- Construct (read) given an IOobject
        featureEdgeMesh(const IOobject&);

        //- Construct from edgeMesh data
        featureEdgeMesh
        (
            const IOobject&,
            const pointField&,
            const edgeList&
        );

        //- Construct as copy
        featureEdgeMesh(const IOobject&, const edgeMesh&);

    //- Destructor
    ~featureEdgeMesh();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "featureEdgeMeshI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
