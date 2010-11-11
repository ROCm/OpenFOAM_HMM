/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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
    Foam::ensightSurfaceWriter

Description

SourceFiles
    ensightSurfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef ensightSurfaceWriter_H
#define ensightSurfaceWriter_H

#include "surfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ensightSurfaceWriter Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class ensightSurfaceWriter
:
    public surfaceWriter<Type>
{
    // Private Member Functions

        static void binShapes
        (
            const pointField& points,
            const faceList& faces,
            DynamicList<label>& tris,
            DynamicList<label>& quads,
            DynamicList<label>& polys
        );

        static void writeGeometry
        (
            const fileName&,
            Ostream&,
            const pointField&,
            const faceList&,
            const DynamicList<label>&,
            const DynamicList<label>&,
            const DynamicList<label>&
        );

        static void writeData(Ostream&, const Field<Type>& values);


public:

    //- Runtime type information
    TypeName("ensight");


    // Constructors

        //- Construct null
        ensightSurfaceWriter();


    //- Destructor
    virtual ~ensightSurfaceWriter();


    // Member Functions

        // Write

        //- Write geometry to file.
        virtual void write
        (
            const fileName& outputDir,
            const fileName& surfaceName,
            const pointField& points,
            const faceList& faces,
            const bool verbose = false
        ) const;


        //- Writes single surface to file.
        virtual void write
        (
            const fileName& outputDir,
            const fileName& surfaceName,
            const pointField& points,
            const faceList& faces,
            const word& fieldName,
            const Field<Type>& values,
            const bool isNodeValues,
            const bool verbose = false
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ensightSurfaceWriter.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
