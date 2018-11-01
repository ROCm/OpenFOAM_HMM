/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2018 OpenCFD Ltd.
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
    Foam::proxySurfaceWriter

Description
    A surfaceWriter that writes the geometry via the MeshedSurfaceProxy, but
    which does not support any fields.

    \heading Output file locations

    The \c rootdir normally corresponds to something like
    \c postProcessing/\<name\>

    \subheading Geometry
    \verbatim
    rootdir
    `-- timeName
        `-- surfaceName.{obj|stl|..}
    \endverbatim

Note
    The formatOptions for proxy are file-type dependent.

SourceFiles
    proxySurfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef proxySurfaceWriter_H
#define proxySurfaceWriter_H

#include "surfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class proxySurfaceWriter Declaration
\*---------------------------------------------------------------------------*/

class proxySurfaceWriter
:
    public surfaceWriter
{
    // Private data

        //- The file extension associated with the proxy
        word fileExtension_;

        //- Format options
        dictionary options_;


public:

    //- Runtime type information
    TypeName("proxy");


    // Constructors

        //- Construct for a given extension
        explicit proxySurfaceWriter(const word& fileExt);

        //- Construct for a given extension, with some output options
        proxySurfaceWriter(const word& fileExt, const dictionary& options);


    //- Destructor
    virtual ~proxySurfaceWriter() = default;


    // Member Functions


        //- True if the surface format supports geometry in a separate file.
        //  False if geometry and field must be in a single file
        virtual bool separateGeometry() const
        {
            return true;
        }


        //- Write single surface geometry to file.
        virtual fileName write
        (
            const fileName& outputDir,      //!< output-dir
            const fileName& surfaceName,    //!< Name of surface
            const meshedSurf& surf,         //!< Surface geometry
            const bool verbose = false      //!< Additional verbosity
        ) const; // override
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
