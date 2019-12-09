/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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
    Foam::ensightPart

Description
    Base class for ensightPartCells and ensightPartFaces

SourceFiles
    ensightPart.C

\*---------------------------------------------------------------------------*/

#ifndef ensightPart_H
#define ensightPart_H

#include "ensightGeoFile.H"
#include "typeInfo.H"
#include "labelList.H"
#include "polyMesh.H"
#include "Field.H"
#include "IOstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ensightPart Declaration
\*---------------------------------------------------------------------------*/

class ensightPart
{
    // Private Data

        //- Part name (or description)
        string name_;


    // Private Member Functions

        //- No copy construct
        ensightPart(const ensightPart&) = delete;

        //- No copy assignment
        void operator=(const ensightPart&) = delete;


protected:

    // Protected Classes

        //- Track the points used by the part and map global to local indices
        struct localPoints
        {
            //- Number of points used
            label nPoints;

            //- Map global to local indices
            labelList list;

            //- Null constructor
            localPoints()
            :
                nPoints(0),
                list()
            {}

            //- Construct for mesh points
            localPoints(const pointField& pts)
            :
                nPoints(0),
                list(pts.size(), -1)
            {}
        };

public:

    //- Runtime type information
    TypeName("ensightPart");


    // Constructors

        //- Construct with description
        explicit ensightPart(const string& description);


    //- Destructor
    virtual ~ensightPart() = default;


    // Access

        //- Part index (0-based)
        virtual label index() const = 0;

        //- Number of elements in this part
        virtual label size() const
        {
            return 0;
        }

        //- Return the part name or description
        const string& name() const
        {
            return name_;
        }

        //- Change the part name or description
        void rename(string value)
        {
            name_ = std::move(value);
        }


    // Output

        //- Write summary information about the object
        virtual void writeSummary(Ostream& os) const = 0;

        //- Print various types of debugging information
        virtual void dumpInfo(Ostream& os) const = 0;

        //- Write geometry
        virtual void write(ensightGeoFile& os) const = 0;

        //- Helper: write geometry with given pointField
        virtual void write(ensightGeoFile& os, const pointField&) const = 0;


    // Housekeeping

        //- Deprecated(2019-12) - use rename() method
        //  \deprecated(2019-12) - use rename() method
        void FOAM_DEPRECATED_FOR(2019-12, "rename() method")
        name(string value)
        {
            name_ = std::move(value);
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


//- IOstream Operator to write geometry
ensightGeoFile& operator<<(ensightGeoFile& os, const ensightPart& part);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
