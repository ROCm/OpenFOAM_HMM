/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::labelToFace

Description
    A topoSetFaceSource to select faces given explicitly provided face labels.

    \heading Dictionary parameters
    \table
        Property    | Description                       | Required  | Default
        normal      | The normal for selecting faces    | yes   |
        cos         | Tolerance angle (range -1, +1)    | yes   |
    \endtable

SourceFiles
    labelToFace.C

\*---------------------------------------------------------------------------*/

#ifndef labelToFace_H
#define labelToFace_H

#include "topoSetFaceSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class labelToFace Declaration
\*---------------------------------------------------------------------------*/

class labelToFace
:
    public topoSetFaceSource
{

    // Private data

        //- Add usage string
        static addToUsageTable usage_;

        //- Cell labels read from dictionary
        labelList labels_;


public:

    //- Runtime type information
    TypeName("labelToFace");

    // Constructors

        //- Construct from components, copying labels
        labelToFace(const polyMesh& mesh, const labelList& labels);

        //- Construct from components, moving labels
        labelToFace(const polyMesh& mesh, labelList&& labels);

        //- Construct from dictionary
        labelToFace(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        labelToFace(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~labelToFace() = default;


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
