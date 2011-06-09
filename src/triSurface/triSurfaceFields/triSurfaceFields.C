/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "triSurfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#ifndef __clang__
template<>
#endif
const word triSurfaceLabelField::typeName("triSurfaceLabelField");

#ifndef __clang__
template<>
#endif
const word triSurfaceScalarField::typeName("triSurfaceScalarField");

#ifndef __clang__
template<>
#endif
const word triSurfaceVectorField::typeName("triSurfaceVectorField");

#ifndef __clang__
template<>
#endif
const word triSurfaceSphericalTensorField::typeName
("triSurfaceSphericalTensorField");

#ifndef __clang__
template<>
#endif
const word triSurfaceSymmTensorField::typeName
("triSurfaceSymmTensorField");

#ifndef __clang__
template<>
#endif
const word triSurfaceTensorField::typeName("triSurfaceTensorField");

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
