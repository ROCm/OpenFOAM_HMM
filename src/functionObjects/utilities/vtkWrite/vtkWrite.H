/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::functionObjects::vtkWrite

Group
    grpUtilitiesFunctionObjects

Description
    This functionObject writes objects registered to the database in VTK format.

    Currently only supports writing of the cell-values of volFields but
    support for other field types, patch fields, Lagrangian data etc. will be
    added.

    Example of function object specification:
    \verbatim
    vtkWrite1
    {
        type            vtkWrite;
        libs            ("libutilityFunctionObjects.so");
        writeControl    writeTime;
        writeInterval   1;
        format          binary;
        legacy          false;
        ...
        fields          (U p);
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: vtkWrite     | yes         |
        fields       | Fields to output        | yes         |
        writeControl | Output control          | recommended | timeStep
        directory    | The output directory name | no        | "VTK"
        format       | ASCII or binary format  | no          | binary
        legacy       | Legacy VTK output       | no          | false
        writeIds     | Write cell ids as field | no          | true
    \endtable

See also
    Foam::functionObjects::ensightWrite
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::timeControl

SourceFiles
    vtkWrite.C
    vtkWriteTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_vtkWrite_H
#define functionObjects_vtkWrite_H

#include "fvMeshFunctionObject.H"
#include "foamVtkInternalWriter.H"
#include "foamVtkSurfaceMeshWriter.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class vtkWrite Declaration
\*---------------------------------------------------------------------------*/

class vtkWrite
:
    public fvMeshFunctionObject
{
    // Private data

        //- VTK output options
        vtk::outputOptions writeOpts_;

        //- Name of fields to process
        wordReList selectFields_;

        //- Output directory name
        fileName dirName_;

        //- Write cell ids field
        bool writeIds_;


    // Private Member Functions

        //- Count number of selected fields for GeoField type.
        //  Only needed for legacy vtk format.
        template<class GeoField>
        label countFields() const;

        //- Write selected fields for GeoField type.
        template<class GeoField>
        label writeFields(vtk::internalWriter& writer, bool verbose=true) const;


        //- Write selected fields for GeoField type.
        template<class GeoField>
        label writeFields
        (
            vtk::surfaceMeshWriter& writer,
            bool verbose=true
        ) const;


        //- Disallow default bitwise copy construct
        vtkWrite(const vtkWrite&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const vtkWrite&) = delete;


public:

    //- Runtime type information
    TypeName("vtkWrite");


    // Constructors

        //- Construct from Time and dictionary
        vtkWrite
        (
            const word& name,
            const Time& t,
            const dictionary& dict
        );


    //- Destructor
    virtual ~vtkWrite() = default;


    // Member Functions

        //- Read the vtkWrite specification
        virtual bool read(const dictionary& dict);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write fields
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "vtkWriteTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
