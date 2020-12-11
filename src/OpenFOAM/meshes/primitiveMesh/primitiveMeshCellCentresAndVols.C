/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Description
    Efficient cell-centre calculation using face-addressing, face-centres and
    face-areas.

\*---------------------------------------------------------------------------*/

#include "primitiveMesh.H"
#include "primitiveMeshTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::primitiveMesh::calcCellCentresAndVols() const
{
    if (debug)
    {
        Pout<< "primitiveMesh::calcCellCentresAndVols() : "
            << "Calculating cell centres and cell volumes"
            << endl;
    }

    // It is an error to attempt to recalculate cellCentres
    // if the pointer is already set
    if (cellCentresPtr_ || cellVolumesPtr_)
    {
        FatalErrorInFunction
            << "Cell centres or cell volumes already calculated"
            << abort(FatalError);
    }

    // set the accumulated cell centre to zero vector
    cellCentresPtr_ = new vectorField(nCells());
    vectorField& cellCtrs = *cellCentresPtr_;

    // Initialise cell volumes to 0
    cellVolumesPtr_ = new scalarField(nCells());
    scalarField& cellVols = *cellVolumesPtr_;

    // Make centres and volumes
    primitiveMeshTools::makeCellCentresAndVols
    (
        *this,
        faceCentres(),
        faceAreas(),
        cellCtrs,
        cellVols
    );

    if (debug)
    {
        Pout<< "primitiveMesh::calcCellCentresAndVols() : "
            << "Finished calculating cell centres and cell volumes"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::vectorField& Foam::primitiveMesh::cellCentres() const
{
    if (!cellCentresPtr_)
    {
        //calcCellCentresAndVols();
        const_cast<primitiveMesh&>(*this).updateGeom();
    }

    return *cellCentresPtr_;
}


const Foam::scalarField& Foam::primitiveMesh::cellVolumes() const
{
    if (!cellVolumesPtr_)
    {
        //calcCellCentresAndVols();
        const_cast<primitiveMesh&>(*this).updateGeom();
    }

    return *cellVolumesPtr_;
}


// ************************************************************************* //
