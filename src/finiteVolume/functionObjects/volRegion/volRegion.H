/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
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
    Foam::functionObjects::volRegion

Group
    grpFieldFunctionObjects

Description
    Volume (cell) region selection class.

    Examples of function object specification:
    \verbatim
    volRegion0
    {
        .
        .
        regionType      cellZone;
        name            c0;
        .
        .
    }

    volRegionAll
    {
        .
        .
        regionType      all;
        .
        .
    }
    \endverbatim

Usage
    \table
        Property     | Description                          | Required | Default
        regionType   | Selection type: all/cellSet/cellZone | no | all
        name         | Name of cellSet/cellZone if required | no |
    \endtable

See also
    Foam::functionObject

SourceFiles
    volRegion.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_volRegion_H
#define functionObjects_volRegion_H

#include "writeFile.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class fvMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class volRegion Declaration
\*---------------------------------------------------------------------------*/

class volRegion
{
    // Private Member Data

        const fvMesh& volMesh_;

        //- The cell ids, from cellSet
        labelList cellIds_;

        // Cache integral properties of the region for writeFileHeader
        label nCells_;

        scalar V_;


public:

    // Public data types

        //- Region type enumeration
        enum regionTypes
        {
            vrtAll,             //!< All cells
            vrtCellSet,         //!< A cellSet
            vrtCellZone         //!< A cellZone
        };

        //- Region type names
        static const Enum<regionTypes> regionTypeNames_;


protected:

    // Protected Data

        //- Region type
        regionTypes regionType_;

        //- Region name (cellSet, cellZone, ...)
        word regionName_;

        //- Region ID (zone ID, ...)
        label regionID_;


    // Protected Member Functions

        //- Output file header information
        void writeFileHeader(const writeFile& wf, Ostream& file) const;


public:

    //- Run-time type information
    TypeName("volRegion");


    // Constructors

        //- Construct from fvMesh and dictionary
        volRegion(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~volRegion() = default;


    // Public Member Functions

        //- Read from dictionary
        bool read(const dictionary& dict);

        //- Return the region type
        inline const regionTypes& regionType() const;

        //- Return the local list of cell IDs
        const labelList& cellIDs() const;

        //- Return the number of cells selected in the region
        label nCells() const;

        //- Return total volume of the region
        scalar V() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volRegionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
