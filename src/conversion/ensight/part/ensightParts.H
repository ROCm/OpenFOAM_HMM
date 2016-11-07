/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::ensightParts

Description
    A collection of several ensightPart elements

SourceFiles
    ensightParts.C

\*---------------------------------------------------------------------------*/

#ifndef ensightParts_H
#define ensightParts_H

#include "SLPtrList.H"
#include "ensightPart.H"
#include "ensightPartFaces.H"
#include "ensightPartCells.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ensightParts Declaration
\*---------------------------------------------------------------------------*/

class ensightParts
:
    public SLPtrList<ensightPart>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        ensightParts(const ensightParts&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const ensightParts&) = delete;


public:

    //- Storage type used
    typedef SLPtrList<ensightPart> StorageType;


    // Constructors

        //- Construct from polyMesh
        ensightParts(const polyMesh&);


    //- Destructor
    ~ensightParts();


    // Member functions

        //- Clear old information and construct anew from polyMesh
        void recalculate(const polyMesh&);

        //- Number of parts
        using StorageType::size;


        //- Write the geometry
        void write(ensightGeoFile&) const;

        //- Write summary information about the objects
        void writeSummary(Ostream&) const;

        //- Print various types of debugging information
        void dumpInfo(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
