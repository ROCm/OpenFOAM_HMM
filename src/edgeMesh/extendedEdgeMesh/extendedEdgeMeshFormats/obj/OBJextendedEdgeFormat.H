/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::fileFormats::OBJextendedEdgeFormat

Description
    Provide a means of reading/writing Alias/Wavefront OBJ format.

    Does not handle negative vertex indices.

SourceFiles
    OBJextendedEdgeFormat.C

\*---------------------------------------------------------------------------*/

#ifndef OBJextendedEdgeFormat_H
#define OBJextendedEdgeFormat_H

#include "extendedEdgeMesh.H"
//#include "IFstream.H"
//#include "Ostream.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                        Class OBJextendedEdgeFormat Declaration
\*---------------------------------------------------------------------------*/

class OBJextendedEdgeFormat
:
    public extendedEdgeMesh
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        OBJextendedEdgeFormat(const OBJextendedEdgeFormat&);

        //- Disallow default bitwise assignment
        void operator=(const OBJextendedEdgeFormat&);


public:

    // Constructors

        //- Construct from file name
        OBJextendedEdgeFormat(const fileName&);


//    // Selectors
//
//        //- Read file and return surface
//        static autoPtr<extendedEdgeMesh> New(const fileName& name)
//        {
//            return autoPtr<extendedEdgeMesh>
//            (
//                new OBJextendedEdgeFormat(name)
//            );
//        }
//

    //- Destructor
    virtual ~OBJextendedEdgeFormat()
    {}


    // Member Functions

        //- Write surface mesh components by proxy
        static void write(const fileName&, const extendedEdgeMesh&);

        //- Read from file
        virtual bool read(const fileName&);

        //- Write object file
        virtual void write(const fileName& name) const
        {
            write(name, *this);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
