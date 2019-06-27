/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2019 OpenCFD Ltd.
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
    Foam::vtk::indirectPatchWriter

Description
    Write indirectPrimitivePatch faces/points (optionally with fields)
    as a vtp file or a legacy vtk file.

    The file output states are managed by the Foam::vtk::fileWriter class.
    FieldData (eg, TimeValue) must appear before any geometry pieces.

Note
    Parallel output is combined into a single Piece without point merging,
    which is similar to using multi-piece data sets, but allows more
    convenient creation as a streaming process.
    In the future, the duplicate points at processor connections
    may be addressed using ghost points.

SourceFiles
    foamVtkIndPatchWriter.C
    foamVtkIndPatchWriterTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkIndPatchWriter_H
#define foamVtkIndPatchWriter_H

#include "foamVtkFileWriter.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
                  Class vtk::indirectPatchWriter Declaration
\*---------------------------------------------------------------------------*/

class indirectPatchWriter
:
    public vtk::fileWriter
{
    // Private Member Data

        //- Reference to the faces as an indirect patch
        const indirectPrimitivePatch& pp_;

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
        //- and begin piece
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
        indirectPatchWriter(const indirectPatchWriter&) = delete;

        //- No copy assignment
        void operator=(const indirectPatchWriter&) = delete;


public:

    // Constructors

        //- Construct from patch (default output INLINE_BASE64)
        indirectPatchWriter
        (
            const indirectPrimitivePatch& pp,
            const vtk::outputOptions opts = vtk::formatType::INLINE_BASE64
        );

        //- Construct from components (default output INLINE_BASE64),
        //  and open the file for writing.
        //  The file name is with/without an extension.
        indirectPatchWriter
        (
            const indirectPrimitivePatch& pp,
            const fileName& file,
            bool parallel = Pstream::parRun()
        );

        //- Construct from components and open the file for writing.
        //  The file name is with/without an extension.
        indirectPatchWriter
        (
            const indirectPrimitivePatch& pp,
            const vtk::outputOptions opts,
            const fileName& file,
            bool parallel = Pstream::parRun()
        );


    //- Destructor
    virtual ~indirectPatchWriter() = default;


    // Member Functions

        //- File extension for current format type.
        using vtk::fileWriter::ext;

        //- File extension for given output type
        inline static word ext(vtk::outputOptions opts)
        {
            return opts.ext(vtk::fileTag::POLY_DATA);
        }

        //- Reference to the indirect patch
        inline const indirectPrimitivePatch& patch() const
        {
            return pp_;
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


    // Write

        //- Write a uniform field of Cell (Face) or Point values
        template<class Type>
        void writeUniform(const word& fieldName, const Type& val);

        //- Write a list of Cell (Face) or Point values
        template<class Type>
        void write(const word& fieldName, const UList<Type>& field);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamVtkIndPatchWriterTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
