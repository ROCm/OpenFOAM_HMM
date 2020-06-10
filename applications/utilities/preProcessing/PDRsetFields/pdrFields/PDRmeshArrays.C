/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 Shell Research Ltd.
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "PDRmeshArrays.H"
#include "PDRblock.H"
#include "polyMesh.H"
#include "Time.H"
#include "IjkField.H"

// Notes
//
// Determines the face and cell numbers of all faces and cells in the
// central rectangular region where CAD_PDR operates. First,
// "points" is read and the coordinates (by which I mean here the
// indices in the x, y and z coordinate arrays) are determined. Then
// "faces" is read and for each the coordinates of the lower- x,y,z
// corner are determioned, also the orientation (X, Y or Z).
// (Orientation in the sense of e.g. + or -x is not noted.) The files
// "owner" and "neighbour" specify the six faces around each cell, so
// from these the coordinates of the cells are determined.
//
// Full checks are made that the mesh in the central region is consistent
// with CAD_PDR's mesh specified by the PDRmeshSpec file.
//
// Eventually, when writing out results, we shall work through the
// full list of cells, writing default values for any cells that are
// not in the central regtion.


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::scalar Foam::PDRmeshArrays::gridPointRelTol = 0.02;


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PDRmeshArrays::classify
(
    const polyMesh& mesh,
    const PDRblock& pdrBlock
)
{
    // Additional copy of i-j-k addressing
    cellDims = pdrBlock.sizes();
    faceDims = (cellDims + labelVector::one);

    const label maxPointId = cmptMax(pdrBlock.sizes())+1;

    Info<< "Mesh" << nl
        << "  nPoints:" << mesh.nPoints()
        << "  nCells:" << mesh.nCells()
        << "  nFaces:" << mesh.nFaces() << nl;

    Info<< "PDRblock" << nl
        << "  minEdgeLen:" << pdrBlock.minEdgeLen() << nl;


    // Bin points into i-j-k locations
    List<labelVector> pointIndex(mesh.nPoints());

    for (label pointi=0; pointi < mesh.nPoints(); ++pointi)
    {
        const point& pt = mesh.points()[pointi];
        pointIndex[pointi] = pdrBlock.gridIndex(pt, gridPointRelTol);
    }

    // Min x,y,z index
    const labelMinMax invertedLimits(maxPointId, -maxPointId);
    Vector<labelMinMax> faceLimits;

    const Vector<direction> faceBits
    (
        boundBox::XDIR,
        boundBox::YDIR,
        boundBox::ZDIR
    );

    faceIndex.resize(mesh.nFaces());
    faceOrient.resize(mesh.nFaces());

    for (label facei=0; facei < mesh.nFaces(); ++facei)
    {
        faceLimits.x() = faceLimits.y() = faceLimits.z() = invertedLimits;

        for (const label pointi : mesh.faces()[facei])
        {
            for (direction cmpt=0; cmpt < labelVector::nComponents; ++cmpt)
            {
                faceLimits[cmpt].add(pointIndex[pointi][cmpt]);
            }
        }

        direction inPlane(0u);

        for (direction cmpt=0; cmpt < labelVector::nComponents; ++cmpt)
        {
            const auto& limits = faceLimits[cmpt];

            if (!limits.valid())
            {
                // This should be impossible
                FatalErrorInFunction
                    << "Unexpected search failure for " << facei << " in "
                    << vector::componentNames[cmpt] << "-direction" << nl
                    << exit(FatalError);
            }

            if (limits.min() < 0)
            {
                FatalErrorInFunction
                    << "Face " << facei << " contains non-grid point in "
                    << vector::componentNames[cmpt] << "-direction" << nl
                    << exit(FatalError);
            }
            else if (limits.min() == limits.max())
            {
                // In plane
                inPlane |= faceBits[cmpt];
            }
            else if (limits.min() + 1 != limits.max())
            {
                FatalErrorInFunction
                    << "Face " << facei
                    << " not in " << vector::componentNames[cmpt] << "-plane" << nl
                    << exit(FatalError);
            }
        }

        switch (inPlane)
        {
            case boundBox::XDIR:
                faceOrient[facei] = vector::X;
                break;

            case boundBox::YDIR:
                faceOrient[facei] = vector::Y;
                break;

            case boundBox::ZDIR:
                faceOrient[facei] = vector::Z;
                break;

            default:
                FatalErrorInFunction
                    << "Face " << facei << " not in an x/y/z plane?" << nl
                    << exit(FatalError);
                break;
        }

        faceIndex[facei] =
            labelVector
            (
                faceLimits.x().min(),
                faceLimits.y().min(),
                faceLimits.z().min()
            );
    }


    // Bin cells into i-j-k locations
    cellIndex = std::move(pointIndex);
    cellIndex = labelVector::uniform(maxPointId);
    cellIndex.resize(mesh.nCells(), labelVector::uniform(maxPointId));

    // Option 1: use PDRblock.findCell() method
    if (true)
    {
        const pointField& cc = mesh.cellCentres();

        for (label celli=0; celli < mesh.nCells(); ++celli)
        {
            cellIndex[celli] = pdrBlock.findCell(cc[celli]);
        }
    }

    // Option 2: walk cell faces and use faceIndex information
    if (false)
    {
        for (label celli=0; celli < mesh.nCells(); ++celli)
        {
            labelVector& cellIdx = cellIndex[celli];

            for (const label facei : mesh.cells()[celli])
            {
                cellIdx.x() = min(cellIdx.x(), faceIndex[facei].x());
                cellIdx.y() = min(cellIdx.y(), faceIndex[facei].y());
                cellIdx.z() = min(cellIdx.z(), faceIndex[facei].z());
            }

            if (cmptMin(cellIdx) < 0)
            {
                cellIdx = labelVector(-1,-1,-1);
            }
        }
    }


    // Verify that all i-j-k cells were found
    {
        // This could be more efficient - but we want to be picky
        IjkField<bool> cellFound(pdrBlock.sizes(), false);

        for (label celli=0; celli < cellIndex.size(); ++celli)
        {
            const labelVector& cellIdx = cellIndex[celli];

            if (cmptMin(cellIdx) >= 0)
            {
                cellFound(cellIdx) = true;
            }
        }

        label firstMissing = cellFound.find(false);

        if (firstMissing >= 0)
        {
            FatalErrorInFunction
                << "No cell found for " << pdrBlock.index(firstMissing)
                << " indexing"
                << exit(FatalError);
        }
    }
}


void Foam::PDRmeshArrays::read
(
    const Time& runTime,
    const PDRblock& pdrBlock
)
{
    #include "createPolyMesh.H"
    classify(mesh, pdrBlock);
}


// ************************************************************************* //
