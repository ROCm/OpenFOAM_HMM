/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "regionModel1D.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
    defineTypeNameAndDebug(regionModel1D, 0);
}
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::regionModels::regionModel1D::constructMeshObjects()
{

    nMagSfPtr_.reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                "nMagSf",
                time().timeName(),
                regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            regionMesh(),
            dimensionedScalar("zero", dimArea, 0.0)
        )
    );
}


void Foam::regionModels::regionModel1D::initialise()
{
    if (debug)
    {
        Pout<< "regionModel1D::initialise()" << endl;
    }

    // Calculate boundaryFaceFaces and boundaryFaceCells

    DynamicList<label> faceIDs;
    DynamicList<label> cellIDs;

    label localPyrolysisFaceI = 0;

    const polyBoundaryMesh& rbm = regionMesh().boundaryMesh();

    forAll(intCoupledPatchIDs_, i)
    {
        const label patchI = intCoupledPatchIDs_[i];
        const polyPatch& ppCoupled = rbm[patchI];
        localPyrolysisFaceI += ppCoupled.size();
    }

    boundaryFaceOppositeFace_.setSize(localPyrolysisFaceI);
    boundaryFaceFaces_.setSize(localPyrolysisFaceI);
    boundaryFaceCells_.setSize(localPyrolysisFaceI);

    localPyrolysisFaceI = 0;

    forAll(intCoupledPatchIDs_, i)
    {
        const label patchI = intCoupledPatchIDs_[i];
        const polyPatch& ppCoupled = rbm[patchI];
        forAll(ppCoupled, localFaceI)
        {
            label faceI = ppCoupled.start() + localFaceI;
            label cellI = -1;
            label nCells = 0;
            do
            {
                label ownCellI = regionMesh().faceOwner()[faceI];
                if (ownCellI != cellI)
                {
                    cellI = ownCellI;
                }
                else
                {
                    cellI = regionMesh().faceNeighbour()[faceI];
                }
                nCells++;
                cellIDs.append(cellI);
                const cell& cFaces = regionMesh().cells()[cellI];
                faceIDs.append(faceI);
                label face0 =
                    cFaces.opposingFaceLabel(faceI, regionMesh().faces());
                faceI = face0;
            } while (regionMesh().isInternalFace(faceI));

            boundaryFaceOppositeFace_[localPyrolysisFaceI] = faceI;
            //faceIDs.remove(); //remove boundary face.

            boundaryFaceFaces_[localPyrolysisFaceI].transfer(faceIDs);
            boundaryFaceCells_[localPyrolysisFaceI].transfer(cellIDs);

            localPyrolysisFaceI++;
            nLayers_ = nCells;
        }
    }
    faceIDs.clear();
    cellIDs.clear();

    surfaceScalarField& nMagSf = nMagSfPtr_();

    localPyrolysisFaceI = 0;
    forAll(intCoupledPatchIDs_, i)
    {
        const label patchI = intCoupledPatchIDs_[i];
        const polyPatch& ppCoupled = rbm[patchI];
        const vectorField& pNormals = ppCoupled.faceNormals();

        nMagSf.boundaryField()[patchI] =
            regionMesh().Sf().boundaryField()[patchI] & pNormals;

        forAll(pNormals, localFaceI)
        {
            const vector n = pNormals[localFaceI];
            const labelList& faces = boundaryFaceFaces_[localPyrolysisFaceI++];

            forAll (faces, faceI) //faceI = 0 is on boundary
            {
                if (faceI > 0)
                {
                    const label faceID = faces[faceI];
                    nMagSf[faceID] = regionMesh().Sf()[faceID] & n;
                }
            }
        }
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::regionModels::regionModel1D::read()
{
    if (regionModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::regionModels::regionModel1D::read(const dictionary& dict)
{
    if (regionModel::read(dict))
    {
        moveMesh_.readIfPresent("moveMesh", coeffs_);

        return true;
    }
    else
    {
        return false;
    }
}


Foam::tmp<Foam::labelField> Foam::regionModels::regionModel1D::moveMesh
(
    const scalarList& deltaV,
    const scalar minDelta
)
{
    tmp<labelField> tcellMoveMap(new labelField(regionMesh().nCells(), 0));
    labelField& cellMoveMap = tcellMoveMap();

    if (!moveMesh_)
    {
        return cellMoveMap;
    }

    pointField oldPoints = regionMesh().points();
    pointField newPoints = oldPoints;

    const polyBoundaryMesh& bm = regionMesh().boundaryMesh();

    label totalFaceId = 0;
    forAll(intCoupledPatchIDs_, localPatchI)
    {
        label patchI = intCoupledPatchIDs_[localPatchI];
        const polyPatch& pp = bm[patchI];

        forAll (pp, patchFaceI)
        {
            const labelList& faces = boundaryFaceFaces_[totalFaceId];
            const labelList& cells = boundaryFaceCells_[totalFaceId];
            const label oFace = boundaryFaceOppositeFace_[totalFaceId];

            const vector n = pp.faceNormals()[patchFaceI];
            const vector sf = pp.faceAreas()[patchFaceI];

            List<point> oldCf(faces.size() + 1, vector::zero);
            List<bool> frozen(faces.size(), false);

            forAll (faces, i)
            {
                oldCf[i] = regionMesh().faceCentres()[faces[i]];
            }

            oldCf[faces.size()] = regionMesh().faceCentres()[oFace];

            forAll (faces, i)
            {
                const label cellI = cells[i];

                if (mag(oldCf[i + 1] - oldCf[i]) < minDelta)
                {
                    frozen[i] = true;
                    cellMoveMap[cellI] = 1;
                }
            }

            vectorField newDelta(cells.size() + 1, vector::zero);

            label j = 0;
            forAllReverse (cells, i)
            {
                const label cellI = cells[i];
                newDelta[j+1] = (deltaV[cellI]/mag(sf))*n + newDelta[j];
                j++;
            }

            forAll (faces, i)
            {
                const label faceI = faces[i];
                const face f = regionMesh().faces()[faceI];

                forAll(f, pti)
                {
                    const label pointI = f[pti];

                    if (!frozen[i])
                    {
                        newPoints[pointI] =
                            oldPoints[pointI] + newDelta[newDelta.size() - 1 - i];
                    }
                }
            }

            totalFaceId ++;
        }
    }
    // Move points
    regionMesh().movePoints(newPoints);

    return tcellMoveMap;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionModels::regionModel1D::regionModel1D
(
    const fvMesh& mesh,
    const word& regionType
)
:
    regionModel(mesh, regionType),
    boundaryFaceFaces_(),
    boundaryFaceCells_(),
    boundaryFaceOppositeFace_(),
    nLayers_(0),
    nMagSfPtr_(NULL),
    moveMesh_(false)
{}


Foam::regionModels::regionModel1D::regionModel1D
(
    const fvMesh& mesh,
    const word& regionType,
    const word& modelName,
    bool readFields
)
:
    regionModel(mesh, regionType, modelName, false),
    boundaryFaceFaces_(regionMesh().nCells()),
    boundaryFaceCells_(regionMesh().nCells()),
    boundaryFaceOppositeFace_(regionMesh().nCells()),
    nLayers_(0),
    nMagSfPtr_(NULL),
    moveMesh_(false)
{
    if (active_)
    {
        constructMeshObjects();
        initialise();
        if (readFields)
        {
            moveMesh_.readIfPresent("moveMesh", coeffs_);
        }
    }
}


Foam::regionModels::regionModel1D::regionModel1D
(
    const fvMesh& mesh,
    const word& regionType,
    const word& modelName,
    const dictionary& dict,
    bool readFields
)
:
    regionModel(mesh, regionType, modelName, dict, readFields),
    boundaryFaceFaces_(regionMesh().nCells()),
    boundaryFaceCells_(regionMesh().nCells()),
    boundaryFaceOppositeFace_(regionMesh().nCells()),
    nLayers_(0),
    nMagSfPtr_(NULL),
    moveMesh_(false)
{
    if (active_)
    {
        constructMeshObjects();
        initialise();
        if (readFields)
        {
            moveMesh_.readIfPresent("moveMesh", coeffs_);
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionModels::regionModel1D::~regionModel1D()
{}


// ************************************************************************* //
