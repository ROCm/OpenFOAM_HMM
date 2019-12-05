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
    Foam::boxToPoint

Description
    A topoSetPointSource to select points based on whether they are inside box.

    \heading Dictionary parameters
    \table
        Property    | Description                           | Required | Default
        boxes       | Multiple bounding boxes               | partly   |
        box         | A single bounding box                 | partly   |
        min         | Minimum point for a single box        | partly   |
        max         | Maximum point for a single box        | partly   |
    \endtable

Note
    Must specify "boxes", "box" or a "min/max" pair
    (compatibility with searchable box) - highest to lowest precedence.

SourceFiles
    boxToPoint.C

\*---------------------------------------------------------------------------*/

#ifndef boxToPoint_H
#define boxToPoint_H

#include "topoSetPointSource.H"
#include "treeBoundBoxList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class boxToPoint Declaration
\*---------------------------------------------------------------------------*/

class boxToPoint
:
    public topoSetPointSource
{

    // Private data

        //- Add usage string
        static addToUsageTable usage_;

        //- Bounding boxes
        treeBoundBoxList bbs_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("boxToPoint");

    // Constructors

        //- Construct from components, copying bounding boxes
        boxToPoint(const polyMesh& mesh, const treeBoundBoxList& bbs);

        //- Construct from components, moving bounding boxes
        boxToPoint(const polyMesh& mesh, treeBoundBoxList&& bbs);

        //- Construct from dictionary
        boxToPoint(const polyMesh& mesh, const dictionary& dict);

        //- Construct from stream - a single box.
        boxToPoint(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~boxToPoint() = default;


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
