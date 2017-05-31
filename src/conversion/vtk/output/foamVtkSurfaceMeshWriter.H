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
    Foam::foamVtkOutput::surfaceMeshWriter

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
#include "indirectPrimitivePatch.H"
#include "foamVtkOutputOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class volPointInterpolation;

namespace foamVtkOutput
{

/*---------------------------------------------------------------------------*\
                      Class surfaceMeshWriter Declaration
\*---------------------------------------------------------------------------*/

class surfaceMeshWriter
{
    const indirectPrimitivePatch& pp_;

    autoPtr<foamVtkOutput::formatter> format_;

    std::ofstream os_;


public:

    // Constructors

        //- Construct from components
        surfaceMeshWriter
        (
            const indirectPrimitivePatch& pp,
            const word& name,
            const fileName& baseName,
            const foamVtkOutput::outputOptions outOpts
        );


    //- Destructor
    ~surfaceMeshWriter();


    // Member Functions

        inline std::ofstream& os()
        {
            return os_;
        }

        inline foamVtkOutput::formatter& format()
        {
            return format_();
        }

        void beginCellData(label nFields);

        void endCellData();


        //- Extract face data
        template<class Type>
        tmp<Field<Type>> getFaceField
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>& sfld
        ) const;

        //- Write surfaceFields
        template<class Type>
        void write
        (
            const UPtrList
            <
                const GeometricField<Type, fvsPatchField, surfaceMesh>
            >& sflds
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace foamVtkOutput
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamVtkSurfaceMeshWriterTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
