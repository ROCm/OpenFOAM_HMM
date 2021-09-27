/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::fileFormats::VTKedgeFormat

Description
    Provide a means of writing VTK legacy format.

SourceFiles
    VTKedgeFormat.C

\*---------------------------------------------------------------------------*/

#ifndef VTKedgeFormat_H
#define VTKedgeFormat_H

#include "edgeMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                     Class VTKedgeFormat Declaration
\*---------------------------------------------------------------------------*/

class VTKedgeFormat
:
    public edgeMesh
{
public:

    // Constructors

        //- Read construct from file name
        explicit VTKedgeFormat(const fileName& filename);


    // Selectors

        //- Read file and return surface
        static autoPtr<edgeMesh> New(const fileName& name)
        {
            return autoPtr<edgeMesh>::NewFrom<VTKedgeFormat>(name);
        }


    //- Destructor
    virtual ~VTKedgeFormat() = default;


    // Member Functions

        //- Write surface mesh components by proxy
        static void write(const fileName& name, const edgeMesh&);

        //- Read from file
        virtual bool read(const fileName& name);

        //- Write object to file
        virtual void write(const fileName& name) const
        {
            write(name, *this);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
