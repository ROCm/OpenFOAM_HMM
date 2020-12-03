/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::isoSurfaceBase

Description
    Low-level components common to various iso-surface algorithms.

SourceFiles
    isoSurfaceBase.C

\*---------------------------------------------------------------------------*/

#ifndef isoSurfaceBase_H
#define isoSurfaceBase_H

#include "isoSurfaceParams.H"
#include "scalar.H"
#include "MeshedSurface.H"
#include "MeshedSurfacesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class isoSurfaceBase Declaration
\*---------------------------------------------------------------------------*/

class isoSurfaceBase
:
    public meshedSurface,
    public isoSurfaceParams
{
protected:

    // Protected typedefs for convenience
    typedef meshedSurface Mesh;

    // Protected Data

        //- Iso value
        const scalar iso_;

        //- For every face, the original cell in mesh
        labelList meshCells_;


public:

    // Constructors

        //- Create for specified algorithm type/filter and iso-value
        isoSurfaceBase
        (
            const scalar iso,
            const isoSurfaceParams& params = isoSurfaceParams()
        );


    // Member Functions

        //- The iso-value associated with the surface
        scalar isoValue() const noexcept
        {
            return iso_;
        }

        //- For each face, the original cell in mesh
        const labelList& meshCells() const noexcept
        {
            return meshCells_;
        }

        //- For each face, the original cell in mesh
        labelList& meshCells() noexcept
        {
            return meshCells_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
