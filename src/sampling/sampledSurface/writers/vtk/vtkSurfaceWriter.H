/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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
    Foam::vtkSurfaceWriter

Description
    A surfaceWriter for VTK legacy format.

SourceFiles
    vtkSurfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef vtkSurfaceWriter_H
#define vtkSurfaceWriter_H

#include "surfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class vtkSurfaceWriter Declaration
\*---------------------------------------------------------------------------*/

class vtkSurfaceWriter
:
    public surfaceWriter
{
    // Private data

        const unsigned int writePrecision_;


    // Private Member Functions

        static void writeGeometry(Ostream&, const meshedSurf&);

        template<class Type>
        static void writeData(Ostream&, const Field<Type>&);


        //- Templated write operation
        template<class Type>
        fileName writeTemplate
        (
            const fileName& outputDir,
            const fileName& surfaceName,
            const meshedSurf& surf,
            const word& fieldName,
            const Field<Type>& values,
            const bool isNodeValues,
            const bool verbose
        ) const;


public:

    //- Runtime type information
    TypeName("vtk");


    // Constructors

        //- Construct null
        vtkSurfaceWriter();

        //- Construct with some output options
        vtkSurfaceWriter(const dictionary& options);


    //- Destructor
    virtual ~vtkSurfaceWriter() = default;


    // Member Functions

        //- Write single surface geometry to file.
        virtual fileName write
        (
            const fileName& outputDir,
            const fileName& surfaceName,
            const meshedSurf& surf,
            const bool verbose = false
        ) const; // override


        //- Write scalarField for a single surface to file.
        //  One value per face or vertex (isNodeValues = true)
        virtual fileName write
        (
            const fileName& outputDir,      // <case>/surface/TIME
            const fileName& surfaceName,    // name of surface
            const meshedSurf& surf,
            const word& fieldName,          // name of field
            const Field<scalar>& values,
            const bool isNodeValues,
            const bool verbose = false
        ) const; // override

        //- Write vectorField for a single surface to file.
        //  One value per face or vertex (isNodeValues = true)
        virtual fileName write
        (
            const fileName& outputDir,      // <case>/surface/TIME
            const fileName& surfaceName,    // name of surface
            const meshedSurf& surf,
            const word& fieldName,          // name of field
            const Field<vector>& values,
            const bool isNodeValues,
            const bool verbose = false
        ) const; // override

        //- Write sphericalTensorField for a single surface to file.
        //  One value per face or vertex (isNodeValues = true)
        virtual fileName write
        (
            const fileName& outputDir,      // <case>/surface/TIME
            const fileName& surfaceName,    // name of surface
            const meshedSurf& surf,
            const word& fieldName,          // name of field
            const Field<sphericalTensor>& values,
            const bool isNodeValues,
            const bool verbose = false
        ) const; // override

        //- Write symmTensorField for a single surface to file.
        //  One value per face or vertex (isNodeValues = true)
        virtual fileName write
        (
            const fileName& outputDir,      // <case>/surface/TIME
            const fileName& surfaceName,    // name of surface
            const meshedSurf& surf,
            const word& fieldName,          // name of field
            const Field<symmTensor>& values,
            const bool isNodeValues,
            const bool verbose = false
        ) const; // override

        //- Write tensorField for a single surface to file.
        //  One value per face or vertex (isNodeValues = true)
        virtual fileName write
        (
            const fileName& outputDir,      // <case>/surface/TIME
            const fileName& surfaceName,    // name of surface
            const meshedSurf& surf,
            const word& fieldName,          // name of field
            const Field<tensor>& values,
            const bool isNodeValues,
            const bool verbose = false
        ) const; // override

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "vtkSurfaceWriterTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
