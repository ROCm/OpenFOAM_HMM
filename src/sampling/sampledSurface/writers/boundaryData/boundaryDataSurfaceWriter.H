/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenCFD Ltd.
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
    Foam::boundaryDataSurfaceWriter

Description
    A surfaceWriter for outputting to a form useable for the
    timeVaryingMapped boundary condition. This reads the data from
    constant/boundaryData/\<patch\> directory

    Typical way of working:
    - use a sampledSurface of type 'patch' (to sample a patch):
    \verbatim
    surfaces
    {
        type            surfaces;
        surfaceFormat   boundaryData;
        fields          ( p );
        surfaces
        (
            outlet
            {
                type            patch;
                patches         (outlet);
                interpolate     false;
            }
        );
    }
    \endverbatim

    - write using this writer.
    - move postProcessing/surfaces/outlet to constant/boundaryData/outlet
      in your destination case.
    - use a timeVaryingMappedFixedValue condition to read and interpolate
      the profile:
        type            timeVaryingMappedFixedValue;
        setAverage      false;  // do not use read average
        offset          0;      // do not apply offset to values

    Note:
    - with 'interpolate false' the data is on the face centres of the
      patch. Take care that a 2D geometry will only have a single row
      of face centres so might not provide a valid triangulation
      (this is what timeVaryingMappedFixedValue uses to do interpolation)
      (Alternatively use timeVaryingMappedFixedValue with mapMethod 'nearest')


SourceFiles
    boundaryDataSurfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryDataSurfaceWriter_H
#define boundaryDataSurfaceWriter_H

#include "surfaceWriter.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class boundaryDataSurfaceWriter Declaration
\*---------------------------------------------------------------------------*/

class boundaryDataSurfaceWriter
:
    public surfaceWriter
{
    // Private Member Functions

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
    TypeName("boundaryData");


    // Constructors

        //- Construct null
        boundaryDataSurfaceWriter();


    //- Destructor
    virtual ~boundaryDataSurfaceWriter();


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
    #include "boundaryDataSurfaceWriterTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
