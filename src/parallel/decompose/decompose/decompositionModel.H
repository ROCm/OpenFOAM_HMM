/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2014-2016 OpenFOAM Foundation
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
    Foam::decompositionModel

Description
    MeshObject wrapper of decompositionMethod

SourceFiles
    decompositionModel.C

\*---------------------------------------------------------------------------*/

#ifndef decompositionModel_H
#define decompositionModel_H

#include "IOdictionary.H"
#include "MeshObject.H"
#include "decompositionMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class mapPolyMesh;
class polyMesh;

/*---------------------------------------------------------------------------*\
                      Class decompositionModel Declaration
\*---------------------------------------------------------------------------*/

class decompositionModel
:
    public MeshObject
    <
        polyMesh,
        UpdateableMeshObject,
        decompositionModel
    >,
    public IOdictionary
{

    // Private Data

        mutable autoPtr<decompositionMethod> decomposerPtr_;


public:

    // Declare name of the class and its debug switch
    ClassName("decompositionModel");

    //- The canonical name ("decomposeParDict") under which the
    //- MeshObject is registered
    static const word canonicalName;


    // Selectors

        //- Read (optionally from absolute path) and register on mesh
        static const decompositionModel& New
        (
            const polyMesh& mesh,
            const fileName& decompDictFile = ""
        );

        //- Read (optionally from supplied dictionary) and register on mesh
        static const decompositionModel& New
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const fileName& decompDictFile = ""
        );


    // Constructors

        //- Construct from typeName or optional path to controlDictionary
        decompositionModel
        (
            const polyMesh& mesh,
            const fileName& decompDictFile = ""
        );


        //- Construct from typeName or optional path to controlDictionary
        decompositionModel
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const fileName& decompDictFile = ""
        );


    // Member Functions

        decompositionMethod& decomposer() const
        {
            if (!decomposerPtr_.valid())
            {
                decomposerPtr_ =
                    decompositionMethod::New
                    (
                        *this,
                        this->mesh().name()  // Name of mesh region
                    );
            }
            return *decomposerPtr_;
        }


      // UpdateableMeshObject

        virtual bool movePoints()
        {
            return false;
        }

        virtual void updateMesh(const mapPolyMesh&)
        {}


    // Housekeeping

        //- Deprecated(2018-08) compatibility method
        //  \deprecated(2018-08) - use IOobject::selectIO directly
        static IOobject selectIO
        (
            const IOobject& io,
            const fileName& altFile,
            const word& ioName = ""
        )
        {
            return IOobject::selectIO(io, altFile, ioName);
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
