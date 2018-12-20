/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::topoSetPointSource

Description
    Base class of a topoSet source for selecting points.

SourceFiles
    topoSetPointSource.C

\*---------------------------------------------------------------------------*/

#ifndef topoSetPointSource_H
#define topoSetPointSource_H

#include "topoSetSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class topoSetPointSource Declaration
\*---------------------------------------------------------------------------*/

class topoSetPointSource
:
    public topoSetSource
{
public:

    //- Runtime type information
    TypeName("topoSetPointSource");

    // Declare run-time constructor selection table

        // For the dictionary constructor
        declareRunTimeSelectionTable
        (
            autoPtr,
            topoSetPointSource,
            word,
            (
                const polyMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );

        // For the Istream constructor
        declareRunTimeSelectionTable
        (
            autoPtr,
            topoSetPointSource,
            istream,
            (
                const polyMesh& mesh,
                Istream& is
            ),
            (mesh, is)
        );


    // Constructors

        //- Construct from components
        explicit topoSetPointSource(const polyMesh& mesh);

        //- Clone (disallowed)
        autoPtr<topoSetPointSource> clone() const
        {
            NotImplemented;
            return nullptr;
        }


    // Selectors

        //- Return a reference to the selected source type
        static autoPtr<topoSetPointSource> New
        (
            const word& sourceType,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Return a reference to the selected source type
        static autoPtr<topoSetPointSource> New
        (
            const word& sourceType,
            const polyMesh& mesh,
            Istream& is
        );


    //- Destructor
    virtual ~topoSetPointSource() = default;


    // Member Functions

        virtual topoSetSource::sourceType setType() const
        {
            return POINTSETSOURCE;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
