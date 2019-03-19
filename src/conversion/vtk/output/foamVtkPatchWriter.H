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
    Foam::vtk::patchWriter

Description
    Write OpenFOAM patches and patch fields in VTP or legacy vtk format.

    The file output states are managed by the Foam::vtk::fileWriter class.
    FieldData (eg, TimeValue) must appear before any geometry pieces.

Note
    Parallel output is combined into a single Piece without point merging,
    which is similar to using multi-piece data sets, but allows more
    convenient creation as a streaming process.
    In the future, the duplicate points at processor connections
    may be addressed using ghost points.

SourceFiles
    foamVtkPatchWriter.C
    foamVtkPatchWriterTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkPatchWriter_H
#define foamVtkPatchWriter_H

#include "foamVtkFileWriter.H"
#include "volFields.H"
#include "pointFields.H"
#include "PrimitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
                      Class vtk::patchWriter Declaration
\*---------------------------------------------------------------------------*/

class patchWriter
:
    public vtk::fileWriter
{
    // Private Member Data

        //- Reference to the OpenFOAM mesh (or subset)
        const fvMesh& mesh_;

        //- The selected patch ids
        labelList patchIDs_;

        //- Use internal field value instead of patch value
        bool useNearCellValue_;

        //- The numer of field points for the current Piece
        label numberOfPoints_;

        //- The numer of field cells (faces) for the current Piece
        label numberOfCells_;

        //- Local number of points
        label nLocalPoints_;

        //- Local number of faces
        label nLocalFaces_;

        //- Local face vertices (connectivity) count. Sum of face sizes.
        label nLocalVerts_;


    // Private Member Functions

        //- Determing sizes (nLocalPoints_, nLocalFaces_, nLocalVerts_),
        //- and begin piece.
        void beginPiece();

        //- Write patch points
        void writePoints();

        //- Write patch faces, legacy format
        //  \param pointOffset processor-local point offset
        void writePolysLegacy(const label pointOffset);

        //- Write patch faces
        //  \param pointOffset processor-local point offset
        void writePolys(const label pointOffset);


        //- No copy construct
        patchWriter(const patchWriter&) = delete;

        //- No copy assignment
        void operator=(const patchWriter&) = delete;


public:

    // Constructors

        //- Construct from components (default format INLINE_BASE64)
        //  \param useNearCellValue to use cell instead of patch values
        patchWriter
        (
            const fvMesh& mesh,
            const labelList& patchIDs,
            const vtk::outputOptions opts = vtk::formatType::INLINE_BASE64,
            const bool useNearCellValue = false
        );

        //- Construct from components (default format INLINE_BASE64),
        //- and open the file for writing.
        //  The file name is with/without an extension.
        patchWriter
        (
            const fvMesh& mesh,
            const labelList& patchIDs,
            const fileName& file,
            bool parallel = Pstream::parRun()
        );

        //- Construct from components (default format INLINE_BASE64),
        //- Construct from components and open the file for writing.
        //  The file name is with/without an extension.
        patchWriter
        (
            const fvMesh& mesh,
            const labelList& patchIDs,
            const vtk::outputOptions opts,
            const fileName& file,
            bool parallel = Pstream::parRun()
        );

        //- Construct from components and open the file for writing.
        //  The file name is with/without an extension.
        //  \param useNearCellValue to use cell instead of patch values
        patchWriter
        (
            const fvMesh& mesh,
            const labelList& patchIDs,
            const vtk::outputOptions opts,
            const bool useNearCellValue,
            const fileName& file,
            bool parallel = Pstream::parRun()
        );


    //- Destructor
    virtual ~patchWriter() = default;


    // Member Functions

        //- File extension for current format type.
        using vtk::fileWriter::ext;

        //- File extension for given output type
        inline static word ext(vtk::outputOptions opts)
        {
            return opts.ext(vtk::fileTag::POLY_DATA);
        }

        //- The patch IDs
        inline const labelList& patchIDs() const
        {
            return patchIDs_;
        }

        //- Write file header (non-collective)
        //  \note Expected calling states: (OPENED).
        virtual bool beginFile(std::string title = "");

        //- Write patch topology
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


        //- Write patch ids as CellData.
        //  Must be called within the CELL_DATA state.
        void writePatchIDs();

        //- Write processor ids as CellData. This is no-op in serial.
        //  Must be called within the CELL_DATA state.
        bool writeProcIDs();

        //- Write processor neighbour ids as CellData. This is no-op in serial.
        //  Must be called within the CELL_DATA state.
        bool writeNeighIDs();


    // Write

        //- Write a uniform field of Cell (Face) or Point values
        template<class Type>
        void writeUniform(const word& fieldName, const Type& val);

        //- Write point field
        template<class Type, template<class> class PatchField>
        void write
        (
            const GeometricField<Type, PatchField, pointMesh>& field
        );

        //- Write volume field
        template<class Type, template<class> class PatchField>
        void write
        (
            const GeometricField<Type, PatchField, volMesh>& field
        );

        //- Write volume field with point interpolation
        template<class Type>
        void write
        (
            const GeometricField<Type, fvPatchField, volMesh>& field,
            const PrimitivePatchInterpolation<primitivePatch>& pInterp
        );

        //- Write volume field with point interpolation
        template<class Type>
        void write
        (
            const GeometricField<Type, fvPatchField, volMesh>& field,
            const PrimitivePatchInterpolation<primitivePatch>* pInterp
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamVtkPatchWriterTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
