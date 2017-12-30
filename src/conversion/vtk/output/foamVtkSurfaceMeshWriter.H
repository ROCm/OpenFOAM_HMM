/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::vtk::surfaceMeshWriter

Description
    Write faces with fields

SourceFiles
    foamVtkSurfaceMeshWriter.C
    foamVtkSurfaceMeshWriterTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkSurfaceMeshWriter_H
#define foamVtkSurfaceMeshWriter_H

#include "pointMesh.H"
#include "OFstream.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "areaFields.H"
#include "indirectPrimitivePatch.H"
#include "foamVtkOutputOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class volPointInterpolation;

namespace vtk
{

/*---------------------------------------------------------------------------*\
                      Class surfaceMeshWriter Declaration
\*---------------------------------------------------------------------------*/

class surfaceMeshWriter
{
    // Private Member Data

        const indirectPrimitivePatch& pp_;

        //- Commonly used query
        const bool legacy_;

        autoPtr<vtk::formatter> format_;

        std::ofstream os_;


    // Private Member Functions

        //- Begin piece
        void beginPiece();

        //- Write patch points
        void writePoints();

        //- Write patch faces
        void writePolysLegacy();

        //- Write patch faces
        void writePolys();

        //- Write mesh topology
        void writeMesh();


        //- Disallow default bitwise copy construct
        surfaceMeshWriter(const surfaceMeshWriter&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const surfaceMeshWriter&) = delete;

public:

    // Constructors

        //- Construct from components
        surfaceMeshWriter
        (
            const indirectPrimitivePatch& pp,
            const word& name,
            const fileName& baseName,
            const vtk::outputOptions outOpts
        );


    //- Destructor
    ~surfaceMeshWriter();


    // Member Functions

        inline std::ofstream& os()
        {
            return os_;
        }

        inline vtk::formatter& format()
        {
            return format_();
        }

        //- Open write for CellData of count fields.
        //  The parameters are only used for the legacy format.
        void beginCellData(label nFields);

        //- Close write for CellData
        void endCellData();

        //- Write file footer
        void writeFooter();


        //- Get face field (internal face or boundary face)
        template<class Type>
        tmp<Field<Type>> getFaceField
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>& sfld
        ) const;


      // Write fields (individually)

        //- Write surface field
        template<class Type>
        void write
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>& field
        );

        //- Write surface field
        template<class Type>
        void write
        (
            const GeometricField<Type, faPatchField, areaMesh>& field
        );


      // Write fields (collectively)

        //- Write surface fields
        template<class Type>
        void write
        (
            const UPtrList
            <
                const GeometricField<Type, fvsPatchField, surfaceMesh>
            >& sflds
        );

        //- Write surface fields
        template<class Type>
        void write
        (
            const UPtrList
            <
                const GeometricField<Type, faPatchField, areaMesh>
            >& sflds
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamVtkSurfaceMeshWriterTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
