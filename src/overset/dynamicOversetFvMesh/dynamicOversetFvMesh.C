/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify i
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

#include "dynamicOversetFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "cellCellStencilObject.H"
#include "zeroGradientFvPatchFields.H"
#include "lduPrimitiveProcessorInterface.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicOversetFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, dynamicOversetFvMesh, IOobject);
}


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicOversetFvMesh::updateAddressing() const
{
    const cellCellStencilObject& overlap = Stencil::New(*this);

    // The (processor-local part of the) stencil determines the local
    // faces to add to the matrix. tbd: parallel
    const labelListList& stencil = overlap.cellStencil();

    // Get the base addressing
    const lduAddressing& baseAddr = dynamicMotionSolverFvMesh::lduAddr();

    // Add to the base addressing
    labelList lowerAddr;
    labelList upperAddr;
    label nExtraFaces;

    const globalIndex globalNumbering(baseAddr.size());
    labelListList localFaceCells;
    labelListList remoteFaceCells;

    labelList globalCellIDs(overlap.cellInterpolationMap().constructSize());
    forAll(baseAddr, cellI)
    {
        globalCellIDs[cellI] = globalNumbering.toGlobal(cellI);
    }
    overlap.cellInterpolationMap().distribute(globalCellIDs);


    reverseFaceMap_ = fvMeshPrimitiveLduAddressing::addAddressing
    (
        baseAddr,
        stencil,
        nExtraFaces,
        lowerAddr,
        upperAddr,
        stencilFaces_,
        globalNumbering,
        globalCellIDs,
        localFaceCells,
        remoteFaceCells
    );

    if (debug)
    {
        Pout<< "dynamicOversetFvMesh::update() : extended addresssing from"
            << " nFaces:" << baseAddr.lowerAddr().size()
            << " to nFaces:" << lowerAddr.size()
            << " nExtraFaces:" << nExtraFaces << endl;
    }

    // Send across wanted cells
    labelListList sendCells;
    Pstream::exchange<labelList, label>(remoteFaceCells, sendCells);

    // Extract relevant remote processors
    labelList nbrProcs(localFaceCells.size());
    {
        label nbrI = 0;
        forAll(localFaceCells, procI)
        {
            if (localFaceCells[procI].size())
            {
                //Pout<< "   from proc:" << procI
                //    << " want its local cells " << remoteFaceCells[procI]
                //    << " to add to my local cells:" << localFaceCells[procI]
                //    << endl;
                nbrProcs[nbrI++] = procI;
            }
        }
        nbrProcs.setSize(nbrI);
    }

    // Construct interfaces
    remoteStencilInterfaces_.setSize(nbrProcs.size());
    forAll(nbrProcs, i)
    {
        label procI = nbrProcs[i];
        remoteStencilInterfaces_.set
        (
            i,
            new lduPrimitiveProcessorInterface
            (
                localFaceCells[procI],
                Pstream::myProcNo(),
                procI,
                tensorField(0),
                Pstream::msgType()
            )
        );
    }


    // Get addressing and interfaces of all interfaces


    List<const labelUList*> patchAddr;
    {
        const fvBoundaryMesh& fvp = boundary();

        patchAddr.setSize(fvp.size() + remoteStencilInterfaces_.size());

        allInterfaces_ = dynamicMotionSolverFvMesh::interfaces();
        allInterfaces_.setSize(patchAddr.size());

        forAll(fvp, patchI)
        {
            patchAddr[patchI] = &fvp[patchI].faceCells();
        }
        forAll(remoteStencilInterfaces_, i)
        {
            label patchI = fvp.size()+i;
            const lduPrimitiveProcessorInterface& pp =
                remoteStencilInterfaces_[i];

            //Pout<< "at patch:" << patchI
            //    << " have procPatch:" << pp.type()
            //    << " from:" << pp.myProcNo()
            //    << " to:" << pp.neighbProcNo()
            //    << " with fc:" << pp.faceCells().size() << endl;

            patchAddr[patchI] = &pp.faceCells();
            allInterfaces_.set(patchI, &pp);
        }
    }
    const lduSchedule ps
    (
        lduPrimitiveMesh::nonBlockingSchedule<processorLduInterface>
        (
            allInterfaces_
        )
    );

    lduPtr_.reset
    (
        new fvMeshPrimitiveLduAddressing
        (
            nCells(),
            lowerAddr.xfer(),
            upperAddr.xfer(),
            patchAddr,
            ps
        )
   );


    // Check
    if (debug)
    {
        const lduAddressing& addr = lduPtr_();  //this->lduAddr();

        Pout<< "Adapted addressing:"
            << " lower:" << addr.lowerAddr().size()
            << " upper:" << addr.upperAddr().size() << endl;

        lduInterfacePtrsList iFaces = this->interfaces();
        // Using lduAddressing::patch
        forAll(patchAddr, patchI)
        {
            Pout<< "    " << patchI << "\tpatchAddr:"
                << addr.patchAddr(patchI).size()
                << endl;
        }

        // Using interfaces
        Pout<< "iFaces:" << iFaces.size() << endl;
        forAll(iFaces, patchI)
        {
            if (iFaces.set(patchI))
            {
                Pout<< "    " << patchI << "\tiFace:" << iFaces[patchI].type()
                        << endl;
            }
        }

        Pout<< "end of printing." << endl;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicOversetFvMesh::dynamicOversetFvMesh(const IOobject& io)
:
    dynamicMotionSolverFvMesh(io),
    active_(false)
{
    // Force loading zoneID field before time gets incremented
    (void)cellCellStencil::zoneID(*this);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicOversetFvMesh::~dynamicOversetFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::lduAddressing& Foam::dynamicOversetFvMesh::lduAddr() const
{
    if (!active_)
    {
        return dynamicMotionSolverFvMesh::lduAddr();
    }
    if (lduPtr_.empty())
    {
        // Build extended addressing
        updateAddressing();
    }
    return lduPtr_();
}


const Foam::fvMeshPrimitiveLduAddressing&
Foam::dynamicOversetFvMesh::primitiveLduAddr() const
{
    if (lduPtr_.empty())
    {
        FatalErrorInFunction
            << "Extended addressing not allocated" << abort(FatalError);
    }

    return lduPtr_();
}


bool Foam::dynamicOversetFvMesh::update()
{
    if (dynamicMotionSolverFvMesh::update())
    {
        // Calculate the local extra faces for the interpolation. Note: could
        // let demand-driven lduAddr() trigger it but just to make sure.
        updateAddressing();

        return true;
    }

    return false;
}


bool Foam::dynamicOversetFvMesh::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool valid
) const
{
    bool ok = dynamicMotionSolverFvMesh::writeObject(fmt, ver, cmp, valid);

    // For postprocessing : write cellTypes and zoneID
    {
        const cellCellStencilObject& overlap = Stencil::New(*this);

        const labelUList& cellTypes = overlap.cellTypes();

        volScalarField volTypes
        (
            IOobject
            (
                "cellTypes",
                this->time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            *this,
            dimensionedScalar("zero", dimless, 0.0),
            zeroGradientFvPatchScalarField::typeName
        );

        forAll(volTypes.internalField(), cellI)
        {
            volTypes[cellI] = cellTypes[cellI];
        }
        volTypes.correctBoundaryConditions();
        volTypes.writeObject(fmt, ver, cmp, valid);
    }
    {
        volScalarField volZoneID
        (
            IOobject
            (
                "zoneID",
                this->time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            *this,
            dimensionedScalar("zero", dimless, 0.0),
            zeroGradientFvPatchScalarField::typeName
        );

        const cellCellStencilObject& overlap = Stencil::New(*this);
        const labelIOList& zoneID = overlap.zoneID();

        forAll(zoneID, cellI)
        {
            volZoneID[cellI] = zoneID[cellI];
        }
        volZoneID.correctBoundaryConditions();
        volZoneID.writeObject(fmt, ver, cmp, valid);
    }
    return ok;
}


// ************************************************************************* //
