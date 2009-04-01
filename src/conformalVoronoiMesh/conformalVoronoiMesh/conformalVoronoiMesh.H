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
    Foam::conformalVoronoiMesh

Description

SourceFiles
    conformalVoronoiMeshI.H
    conformalVoronoiMesh.C
    conformalVoronoiMeshIO.C

\*---------------------------------------------------------------------------*/

#ifndef conformalVoronoiMesh_H
#define conformalVoronoiMesh_H

#include ".H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class conformalVoronoiMesh;
Istream& operator>>(Istream&, conformalVoronoiMesh&);
Ostream& operator<<(Ostream&, const conformalVoronoiMesh&);


/*---------------------------------------------------------------------------*\
                          Class conformalVoronoiMesh Declaration
\*---------------------------------------------------------------------------*/

class conformalVoronoiMesh
:
    public baseClassName
{
    // Private data

        //- Description of data_
        dataType data_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        conformalVoronoiMesh(const conformalVoronoiMesh&);

        //- Disallow default bitwise assignment
        void operator=(const conformalVoronoiMesh&);


public:

    // Static data members

        //- Static data staticData
        static const dataType staticData;


    // Constructors

        //- Construct null
        conformalVoronoiMesh();

        //- Construct from components
        conformalVoronoiMesh(const dataType& data);

        //- Construct from Istream
        conformalVoronoiMesh(Istream&);

        //- Construct as copy
        conformalVoronoiMesh(const conformalVoronoiMesh&);


    // Selectors

        //- Select null constructed
        static autoPtr<conformalVoronoiMesh> New();


    //- Destructor
    ~conformalVoronoiMesh();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

        void operator=(const conformalVoronoiMesh&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, conformalVoronoiMesh&);
        friend Ostream& operator<<(Ostream&, const conformalVoronoiMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "conformalVoronoiMeshI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
