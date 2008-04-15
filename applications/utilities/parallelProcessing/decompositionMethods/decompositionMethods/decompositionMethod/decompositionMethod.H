/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    Foam::decompositionMethod

Description
    Abstract base class for decomposition

SourceFiles
    decompositionMethod.C

\*---------------------------------------------------------------------------*/

#ifndef decompositionMethod_H
#define decompositionMethod_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "dictionary.H"
#include "polyMesh.H"
#include "pointField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class decompositionMethod Declaration
\*---------------------------------------------------------------------------*/

class decompositionMethod
{

protected:

    // Protected data

        const dictionary& decompositionDict_;
        label nProcessors_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        decompositionMethod(const decompositionMethod&);
        void operator=(const decompositionMethod&);


public:

    //- Runtime type information
    TypeName("decompositionMethod");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            decompositionMethod,
            dictionary,
            (
                const dictionary& decompositionDict
            ),
            (decompositionDict)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            decompositionMethod,
            dictionaryMesh,
            (
                const dictionary& decompositionDict,
                const polyMesh& mesh
            ),
            (decompositionDict, mesh)
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<decompositionMethod> New
        (
            const dictionary& decompositionDict
        );

        //- Return a reference to the selected turbulence model
        static autoPtr<decompositionMethod> New
        (
            const dictionary& decompositionDict,
            const polyMesh& mesh
        );


    // Constructors

        //- Construct given the decomposition dictionary
        decompositionMethod(const dictionary& decompositionDict)
        :
            decompositionDict_(decompositionDict),
            nProcessors_
            (
                readLabel(decompositionDict.lookup("numberOfSubdomains"))
            )
        {}


    // Destructor

        virtual ~decompositionMethod()
        {}


    // Member Functions

        //- Is method parallel aware (i.e. does it synchronize domains across
        //  proc boundaries)
        virtual bool parallelAware() const = 0;

        //- Return for every coordinate the wanted processor number
        virtual labelList decompose(const pointField&) = 0;

        //- Return for every coordinate the wanted processor number. Gets
        //  passed agglomeration map (from fine to coarse cells) and coarse cell
        //  location. Can be overridden by decomposers that provide this
        //  functionality natively.
        virtual labelList decompose
        (
            const labelList& agglom,
            const pointField&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
