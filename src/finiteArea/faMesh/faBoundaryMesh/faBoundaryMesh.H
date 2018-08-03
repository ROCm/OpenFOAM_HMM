/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2016-2017 Wikki Ltd
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
    Foam::faBoundaryMesh

Description
    Finite area boundary mesh

SourceFiles
    faBoundaryMesh.C

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#ifndef faBoundaryMesh_H
#define faBoundaryMesh_H

#include "faPatchList.H"
#include "lduInterfacePtrsList.H"
#include "wordList.H"
#include "pointField.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward declarations
class faMesh;
class faBoundaryMesh;
Ostream& operator<<(Ostream&, const faBoundaryMesh&);

/*---------------------------------------------------------------------------*\
                      Class faBoundaryMesh Declaration
\*---------------------------------------------------------------------------*/

class faBoundaryMesh
:
    public faPatchList,
    public regIOobject
{
    // Private Data

        //- Reference to mesh
        const faMesh& mesh_;

        //- No copy construct
        faBoundaryMesh(const faBoundaryMesh&) = delete;

        //- No copy assignment
        void operator=(const faBoundaryMesh&) = delete;


public:

    //- Runtime type information
    TypeName("faBoundaryMesh");


    // Constructors

        //- Construct from dictionary
        faBoundaryMesh
        (
            const IOobject& io,
            const faMesh& fam
        );

        //- Construct given size
        faBoundaryMesh
        (
            const IOobject& io,
            const faMesh& fam,
            const label size
        );


    //- Destructor
    ~faBoundaryMesh() = default;


    // Member Functions

        // Access

            //- Calculate the geometry for the patches
            //  (transformation tensors etc.)
            void calcGeometry();

            //- Return the mesh reference
            const faMesh& mesh() const;

            //- Return a list of pointers for each patch
            //  with only those pointing to interfaces being set
            lduInterfacePtrsList interfaces() const;

            //- Return a list of patch names
            wordList names() const;

            //- Return a list of patch types
            wordList types() const;

            //- Return patch indices for all matches.
            //  Matching patchGroups currently not supported
            labelList indices
            (
                const keyType& key,
                const bool useGroups = false
            ) const;

            //- Return patch index for the first match, return -1 if not found
            label findIndex(const keyType& key) const;

            //- Find patch index given a name, return -1 if not found
            label findPatchID(const word& patchName) const;

            //- Return patch index for a given edge label
            label whichPatch(const label edgeIndex) const;

            //- Check boundary definition
            bool checkDefinition(const bool report = false) const;


        // Edit

            //- Correct faBoundaryMesh after moving points
            void movePoints(const pointField&);

            //- Correct faBoundaryMesh after topology update
            void updateMesh();

            //- writeData member function required by regIOobject
            bool writeData(Ostream&) const;


    // Ostream operator

        friend Ostream& operator<<(Ostream&, const faBoundaryMesh&);


    // Housekeeping

        //- Identical to the indices() method (AUG-2018)
        labelList findIndices
        (
            const keyType& key,
            const bool useGroups = false
        ) const
        {
            return indices(key, useGroups);
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
