/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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
    Foam::foamVtkOutput::internalWriter

Description
    Write fields (internal).

SourceFiles
    foamVtkInternalWriter.C
    foamVtkInternalWriterTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkInternalWriter_H
#define foamVtkInternalWriter_H

#include "OFstream.H"
#include "volFields.H"
#include "pointFields.H"

#include "foamVtkCells.H"
#include "foamVtkOutputFields.H"
#include "foamVtkOutputOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class volPointInterpolation;

namespace foamVtkOutput
{

/*---------------------------------------------------------------------------*\
                       Class internalWriter Declaration
\*---------------------------------------------------------------------------*/

class internalWriter
{
    // Private Member Data

        //- Reference to the OpenFOAM mesh (or subset)
        const fvMesh& mesh_;

        //- Commonly used query
        const bool legacy_;

        autoPtr<foamVtkOutput::formatter> format_;

        //- The volume cells (internalMesh)
        const foamVtkCells& vtkCells_;

        std::ofstream os_;


    // Private Member Functions

        //- Begin piece
        void beginPiece();

        //- Write mesh points
        void writePoints();

        //- Write mesh cells
        void writeCellsLegacy();

        //- Write mesh cells
        void writeCells();

        //- Write mesh topology
        void writeMesh();


        //- Disallow default bitwise copy construct
        internalWriter(const internalWriter&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const internalWriter&) = delete;


public:

    // Constructors

        //- Construct from components
        internalWriter
        (
            const fvMesh& mesh,
            const foamVtkCells& cells,
            const fileName& baseName,
            const foamVtkOutput::outputOptions outOpts
        );


    //- Destructor
    ~internalWriter();


    // Member Functions

        inline std::ofstream& os()
        {
            return os_;
        }

        inline foamVtkOutput::formatter& format()
        {
            return format_();
        }


        //- Open write for CellData of count fields.
        //  The parameters are only used for the legacy format.
        void beginCellData(label nFields);

        //- Close write for CellData
        void endCellData();

        //- Open write for PointData of count fields
        //  The parameters are only used for the legacy format.
        void beginPointData(label nFields);

        //- Close write for PointData
        void endPointData();


        //- Write cellIDs
        void writeCellIDs();

        //- Write file footer
        void writeFooter();


        //- Write internal fields
        template<class Type>
        void write
        (
            const UPtrList
            <
                const DimensionedField<Type, volMesh>
            >& flds
        );

        //- Write volFields
        template<class Type, template<class> class PatchField>
        void write
        (
            const UPtrList
            <
                const GeometricField<Type, PatchField, volMesh>
            >& flds
        );

        //- Write pointFields on all mesh points.
        //  Interpolate to cell centre for decomposed cell centres.
        template<class Type, template<class> class PatchField>
        void write
        (
            const UPtrList
            <
                const GeometricField<Type, PatchField, pointMesh>
            >& flds
        );


        //- Interpolated internal fields
        template<class Type>
        void write
        (
            const volPointInterpolation& pInterp,
            const UPtrList
            <
                const DimensionedField<Type, volMesh>
            >& flds
        );

        //- Interpolated volFields
        template<class Type, template<class> class PatchField>
        void write
        (
            const volPointInterpolation& pInterp,
            const UPtrList
            <
                const GeometricField<Type, PatchField, volMesh>
            >& flds
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace foamVtkOutput
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamVtkInternalWriterTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
