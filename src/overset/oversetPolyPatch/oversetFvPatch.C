/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "oversetFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "transform.H"
#include "cellCellStencilObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(oversetFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, oversetFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::labelField> Foam::oversetFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return patchInternalField(internalData);
}


Foam::tmp<Foam::labelField> Foam::oversetFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& restrictMap
) const
{
    // Store the restrictMap. This routine gets used for
    // - GAMGAgglomeration      : this is the use we want to intercept.
    // - GAMGProcAgglomeration  : to find out the cell number on the other side
    // - MGridGenGAMGAgglomeration: same
    if (master())
    {
        restrictMap_ = restrictMap;
    }

    return patchInternalField(restrictMap);
}


const Foam::labelListList& Foam::oversetFvPatch::stencil() const
{
    const cellCellStencilObject& overlap = Stencil::New(boundaryMesh().mesh());

    return overlap.cellStencil();
}


const Foam::mapDistribute& Foam::oversetFvPatch::cellInterpolationMap() const
{
    const cellCellStencilObject& overlap = Stencil::New(boundaryMesh().mesh());
    return overlap.cellInterpolationMap();
}


const Foam::List<Foam::scalarList>&
Foam::oversetFvPatch::cellInterpolationWeights() const
{
    const cellCellStencilObject& overlap = Stencil::New(boundaryMesh().mesh());
    return overlap.cellInterpolationWeights();
}


const Foam::scalarField& Foam::oversetFvPatch::normalisation() const
{
    return boundaryMesh().mesh().V().field();
}


const Foam::labelList& Foam::oversetFvPatch::interpolationCells() const
{
    const cellCellStencilObject& overlap = Stencil::New(boundaryMesh().mesh());
    return overlap.interpolationCells();
}


const Foam::scalarList& Foam::oversetFvPatch::cellInterpolationWeight() const
{
    const cellCellStencilObject& overlap = Stencil::New(boundaryMesh().mesh());
    return overlap.cellInterpolationWeight();
}


// ************************************************************************* //
