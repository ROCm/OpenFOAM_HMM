/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::ensightOutput

Description
    A collection of functions for writing ensight file content in parallel.

SourceFiles
    ensightOutput.C
    ensightOutputTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef ensightOutput_H
#define ensightOutput_H

#include "ensightFile.H"
#include "ensightMesh.H"
#include "autoPtr.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ensightOutput Declaration
\*---------------------------------------------------------------------------*/

class ensightOutput
{
    // Private Methods

    template<class Type>
    static void writeField
    (
        const char* key,
        const Field<Type>& fld,
        ensightFile& os
    );

    template<class Type>
    static bool writePatchField
    (
        const Field<Type>& pf,
        const ensightFaces& ensFaces,
        ensightFile& os
    );

    template<class Type>
    static bool writeVolField
    (
        const Field<Type>& vf,
        const ensightCells& ensCells,
        ensightFile& os,
        const bool deprecatedOrder = false
    );


    //- Write volume field component-wise
    template<class Type>
    static bool writeField
    (
        const GeometricField<Type, fvPatchField, volMesh>& vf,
        const ensightMesh& ensMesh,
        ensightFile& os
    );

    //- Write point field component-wise
    template<class Type>
    static bool ensightPointField
    (
        const GeometricField<Type, pointPatchField, pointMesh>& pf,
        const ensightMesh& ensMesh,
        ensightFile& os
    );


    //- Disallow null constructor
    ensightOutput() = delete;

public:

    // Public Methods


    //- Write volume field component-wise
    template<class Type>
    static bool writeField
    (
        const GeometricField<Type, fvPatchField, volMesh>&,
        const ensightMesh& ensMesh,
        ensightFile& os,
        const bool nodeValues
    );


    //- Write volume field component-wise
    template<class Type>
    static inline bool writeField
    (
        const GeometricField<Type, fvPatchField, volMesh>& vf,
        const ensightMesh& ensMesh,
        autoPtr<ensightFile>& output,
        const bool nodeValues = false
    )
    {
        return writeField(vf, ensMesh, output.rawRef(), nodeValues);
    }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ensightOutputTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
