/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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
    Foam::DirectInteractionList

Description

SourceFiles
    DirectInteractionListI.H
    DirectInteractionList.C

\*---------------------------------------------------------------------------*/

#ifndef DirectInteractionList_H
#define DirectInteractionList_H

#include "polyMesh.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class ParticleType>
class InteractionLists;

/*---------------------------------------------------------------------------*\
                           Class DirectInteractionList Declaration
\*---------------------------------------------------------------------------*/

template<class ParticleType>
class DirectInteractionList
:
    public labelListList
{
    // Private data

        //- A reference to the InteractionLists object holding this object
        const InteractionLists<ParticleType>& il_;

        //- Data specifying which patch faces of type "wall" are in range
        //  of each cell
        labelListList wallFaces_;


    // Private Member Functions

        void buildDirectInteractionList
        (
            bool pointPointListBuild
        );

        void buildWallFaces();

        //- Disallow default bitwise copy construct
        DirectInteractionList(const DirectInteractionList&);

        //- Disallow default bitwise assignment
        void operator=(const DirectInteractionList&);


public:

    // Constructors

        //- Construct lists by searching the mesh
        DirectInteractionList
        (
            const InteractionLists<ParticleType>& il,
            bool pointPointListBuild
        );

        //- Construct from file
        DirectInteractionList
        (
            const InteractionLists<ParticleType>& il
        );


    // Destructor

        ~DirectInteractionList();


    // Member Functions

        // Access

        inline const labelListList& wallFaces() const;

        inline const InteractionLists<ParticleType>& il() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DirectInteractionListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DirectInteractionList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
