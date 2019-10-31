/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::faceToFace

Description
    A topoSetFaceSource to select faces based on usage in another faceSet.

    \heading Dictionary parameters
    \table
        Property    | Description                       | Required  | Default
        sets        | The input faceSet names           | possibly |
        set         | The input faceSet name            | possibly |
    \endtable

Note
    Must specify "sets" or "set" (highest to lowest precedence).

SourceFiles
    faceToFace.C

\*---------------------------------------------------------------------------*/

#ifndef faceToFace_H
#define faceToFace_H

#include "topoSetFaceSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class faceToFace Declaration
\*---------------------------------------------------------------------------*/

class faceToFace
:
    public topoSetFaceSource
{
    // Private data

        //- Add usage string
        static addToUsageTable usage_;

        //- Names of sets to use
        wordList names_;


public:

    //- Runtime type information
    TypeName("faceToFace");

    // Constructors

        //- Construct from components
        faceToFace(const polyMesh& mesh, const word& setName);

        //- Construct from dictionary
        faceToFace(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        faceToFace(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~faceToFace() = default;


    // Member Functions

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet& set
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
