/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::ensightPartCells

Description
    An implementation of ensightPart to hold volume mesh cells.

SourceFiles
    ensightPartCells.C

\*---------------------------------------------------------------------------*/

#ifndef ensightPartCells_H
#define ensightPartCells_H

#include "ensightPart.H"
#include "ensightCells.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ensightPartCells Declaration
\*---------------------------------------------------------------------------*/

class ensightPartCells
:
    public ensightCells,
    public ensightPart
{
    // Private data

        //- Mesh referenced
        const polyMesh& mesh_;


    // Private Member Functions

        //- Track points used
        localPoints calcLocalPoints() const;

        //- Element connectivity
        void writeConnectivity
        (
            ensightGeoFile&,
            const word& key,
            const labelUList& idList,
            const labelUList& pointMap
        ) const;


        //- Disallow default bitwise copy construct
        ensightPartCells(const ensightPartCells&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const ensightPartCells&) = delete;


public:

    //- Runtime type information
    TypeName("ensightCells");

    // Constructors

        //- Construct from polyMesh without zones
        ensightPartCells
        (
            label partIndex,
            const polyMesh&
        );

        //- Construct from polyMesh and list of (non-zoned) cells
        ensightPartCells
        (
            label partIndex,
            const polyMesh&,
            const labelUList&
        );

        //- Construct from polyMesh and cellZone
        ensightPartCells
        (
            label partIndex,
            const polyMesh&,
            const cellZone&
        );


    //- Destructor
    virtual ~ensightPartCells();


    // Member Functions

    // Access

        //- Part index (0-based)
        virtual label index() const
        {
            return ensightCells::index();
        }

        //- Number of elements in this part
        virtual label size() const
        {
            return ensightCells::size();
        }


    // Output

        //- Write geometry
        virtual void write(ensightGeoFile&) const;

        //- Helper: write geometry given the pointField
        virtual void write(ensightGeoFile&, const pointField&) const;


        //- Write summary information about the object
        virtual void writeSummary(Ostream&) const;

        //- Print various types of debugging information
        virtual void dumpInfo(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
