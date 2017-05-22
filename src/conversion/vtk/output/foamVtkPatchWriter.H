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
    Foam::foamVtkOutput::patchWriter

Description
    Write patch fields

SourceFiles
    foamVtkPatchWriter.C
    foamVtkPatchWriterTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkPatchWriter_H
#define foamVtkPatchWriter_H

#include "pointMesh.H"
#include "OFstream.H"
#include "volFields.H"
#include "pointFields.H"
#include "indirectPrimitivePatch.H"
#include "PrimitivePatchInterpolation.H"
#include "foamVtkOutput.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class volPointInterpolation;

namespace foamVtkOutput
{

/*---------------------------------------------------------------------------*\
                           Class patchWriter Declaration
\*---------------------------------------------------------------------------*/

class patchWriter
{
    //- Reference to the OpenFOAM mesh (or subset)
    const fvMesh& mesh_;

    autoPtr<foamVtkOutput::formatter> format_;

    const bool nearCellValue_;

    const labelList patchIDs_;

    std::ofstream os_;

    label nPoints_;

    label nFaces_;

public:

    // Constructors

        //- Construct from components
        patchWriter
        (
            const fvMesh& mesh,
            const bool binary,
            const bool nearCellValue,
            const fileName& fName,
            const labelList& patchIDs
        );


    // Member Functions

        inline std::ofstream& os()
        {
            return os_;
        }

        inline foamVtkOutput::formatter& format()
        {
            return format_();
        }

        inline label nPoints() const
        {
            return nPoints_;
        }

        inline label nFaces() const
        {
            return nFaces_;
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
        void writePatchIDs();

        //- Write volFields
        template<class Type, template<class> class PatchField>
        void write
        (
            const UPtrList
            <
                const GeometricField<Type, PatchField, volMesh>
            >& flds
        );

        //- Write pointFields
        template<class Type, template<class> class PatchField>
        void write
        (
            const UPtrList
            <
                const GeometricField<Type, PatchField, pointMesh>
            >&
        );

        //- Interpolated volFields
        template<class Type>
        void write
        (
            const PrimitivePatchInterpolation<primitivePatch>& pInterp,
            const UPtrList
            <
                const GeometricField<Type, fvPatchField, volMesh>
            >& flds
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace foamVtkOutput
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamVtkPatchWriterTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
