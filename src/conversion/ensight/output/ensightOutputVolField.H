/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

InNamespace
    Foam::ensightOutput

Description
    A collection of functions for writing volField content in ensight format.

\*---------------------------------------------------------------------------*/

#ifndef ensightOutputVolField_H
#define ensightOutputVolField_H

#include "ensightOutput.H"
#include "ensightFaces.H"
#include "ensightCells.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class ensightMesh;

namespace ensightOutput
{

/*---------------------------------------------------------------------------*\
                            Namespace ensightOutput
\*---------------------------------------------------------------------------*/

//- Write volume field component-wise
template<class Type>
bool writeVolField
(
    ensightFile& os,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const ensightMesh& ensMesh
);


//- Write volume field component-wise, optionally forcing interpolation
template<class Type>
bool writeVolField
(
    ensightFile& os,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const ensightMesh& ensMesh,
    const bool nodeValues
);

//- Write point field component-wise
template<class Type>
bool writePointField
(
    ensightFile& os,
    const GeometricField<Type, pointPatchField, pointMesh>& pf,
    const ensightMesh& ensMesh
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ensightOutput
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ensightOutputVolFieldTemplates.C"
#endif

#endif

// ************************************************************************* //
