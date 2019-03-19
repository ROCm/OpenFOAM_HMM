/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2019 OpenCFD Ltd.
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
    Foam::vtk::internalWriter

Description
    Write an OpenFOAM volume (internal) geometry and internal fields
    as a vtu file or a legacy vtk file.

    The file output states are managed by the Foam::vtk::fileWriter class.
    FieldData (eg, TimeValue) must appear before any geometry pieces.

Note
    Parallel output is combined into a single Piece without point merging,
    which is similar to using multi-piece data sets, but allows more
    convenient creation as a streaming process.
    In the future, the duplicate points at processor connections
    may be addressed using ghost points.

SourceFiles
    foamVtkInternalWriter.C
    foamVtkInternalWriterTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkInternalWriter_H
#define foamVtkInternalWriter_H

#include "foamVtkFileWriter.H"
#include "foamVtuCells.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class volPointInterpolation;

namespace vtk
{

/*---------------------------------------------------------------------------*\
                     Class vtk::internalWriter Declaration
\*---------------------------------------------------------------------------*/

class internalWriter
:
    public vtk::fileWriter
{
    // Private Member Data

        //- Reference to the OpenFOAM mesh (or subset)
        const fvMesh& mesh_;

        //- The volume cells (internalMesh)
        const vtuCells& vtuCells_;

        //- The numer of field points for the current Piece
        label numberOfPoints_;

        //- The numer of field cells for the current Piece
        label numberOfCells_;


    // Private Member Functions

        //- Begin piece
        void beginPiece();

        //- Write mesh points
        void writePoints();

        //- Write cells (connectivity and type), legacy format
        //  \param pointOffset processor-local point offset
        void writeCellsLegacy(const label pointOffset);

        //- Write cells connectivity
        //  \param pointOffset processor-local point offset
        void writeCellsConnectivity(const label pointOffset);

        //- Write cells face streams
        //  \param pointOffset processor-local point offset
        void writeCellsFaces(const label pointOffset);


        //- No copy construct
        internalWriter(const internalWriter&) = delete;

        //- No copy assignment
        void operator=(const internalWriter&) = delete;


public:

    //- Debug information
    static int debug;

//- File extension for given output type
//  constexpr static vtk::fileTag content = vtk::fileTag::UNSTRUCTURED_GRID;


    // Constructors

        //- Construct from components (default format INLINE_BASE64),
        internalWriter
        (
            const fvMesh& mesh,
            const vtk::vtuCells& cells,
            const vtk::outputOptions opts = vtk::formatType::INLINE_BASE64
        );

        //- Construct from components (default format INLINE_BASE64),
        //  The file name is with/without an extension.
        internalWriter
        (
            const fvMesh& mesh,
            const vtk::vtuCells& cells,
            const fileName& file,
            bool parallel = Pstream::parRun()
        );

        //- Construct from components (default format INLINE_BASE64),
        //  The file name is with/without an extension.
        internalWriter
        (
            const fvMesh& mesh,
            const vtk::vtuCells& cells,
            const vtk::outputOptions opts,
            const fileName& file,
            bool parallel = Pstream::parRun()
        );


    //- Destructor
    virtual ~internalWriter() = default;


    // Member Functions

        //- File extension for current format type.
        using vtk::fileWriter::ext;

        //- File extension for given output type
        inline static word ext(vtk::outputOptions opts)
        {
            return opts.ext(vtk::fileTag::UNSTRUCTURED_GRID);
        }


        //- Write file header (non-collective)
        //  \note Expected calling states: (OPENED)
        virtual bool beginFile(std::string title = "");

        //- Write mesh topology.
        //  Also writes the file header if not previously written.
        //  \note Must be called prior to writing CellData or PointData
        virtual bool writeGeometry();

        //- Begin CellData output section for specified number of fields.
        //  Must be called prior to writing any cell data fields.
        //  \param nFields is for legacy format only.
        //      When nFields=0, this a no-op for legacy format.
        //  \note Expected calling states: (PIECE | POINT_DATA).
        //
        //  \return True if the state changed
        virtual bool beginCellData(label nFields = 0);

        //- Begin PointData for specified number of fields.
        //  Must be called prior to writing any point data fields.
        //  \param nFields is for legacy format only.
        //      When nFields=0, this a no-op for legacy format.
        //  \note Expected calling states: (PIECE | CELL_DATA).
        //
        //  \return True if the state changed
        virtual bool beginPointData(label nFields = 0);


        //- Write cell ids as CellData.
        //  Must be called within the CELL_DATA state.
        void writeCellIDs();

        //- Write processor ids as CellData. This is no-op in serial.
        //  Must be called within the CELL_DATA state.
        bool writeProcIDs();

        //- Write point ids as PointData.
        //  Must be called within the POINT_DATA state.
        void writePointIDs();


    // Write

        //- Write a uniform field of Cell or Point values
        template<class Type>
        void writeUniform(const word& fieldName, const Type& val);

        //- Write point field
        //  Interpolate to originating cell centre for decomposed cells.
        template<class Type, template<class> class PatchField>
        void write
        (
            const GeometricField<Type, PatchField, pointMesh>& field
        );

        //- Write the internal field
        template<class Type>
        void write
        (
            const DimensionedField<Type, volMesh>& field
        );

        //- Write the volume field (internal part)
        template<class Type, template<class> class PatchField>
        void write
        (
            const GeometricField<Type, PatchField, volMesh>& field
        );

        //- Write internal field with point interpolation
        template<class Type>
        void write
        (
            const DimensionedField<Type, volMesh>& vfield,
            const volPointInterpolation& pInterp
        );

        //- Write volume field with point interpolation
        template<class Type>
        void write
        (
            const GeometricField<Type, fvPatchField, volMesh>& vfield,
            const volPointInterpolation& pInterp
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamVtkInternalWriterTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
