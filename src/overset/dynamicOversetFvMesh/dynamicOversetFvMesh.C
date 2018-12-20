/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2018 OpenCFD Ltd.
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
        Pout<< "dynamicOversetFvMesh::update() : extended addressing from"
            << " nFaces:" << baseAddr.lowerAddr().size()
            << " to nFaces:" << lowerAddr.size()
            << " nExtraFaces:" << nExtraFaces << endl;
    }

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
            std::move(lowerAddr),
            std::move(upperAddr),
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
    // Load stencil (but do not update)
    (void)Stencil::New(*this, false);
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
    return *lduPtr_;
}


const Foam::fvMeshPrimitiveLduAddressing&
Foam::dynamicOversetFvMesh::primitiveLduAddr() const
{
    if (lduPtr_.empty())
    {
        FatalErrorInFunction
            << "Extended addressing not allocated" << abort(FatalError);
    }

    return *lduPtr_;
}


bool Foam::dynamicOversetFvMesh::update()
{
    if (dynamicMotionSolverFvMesh::update())
    {
        // Calculate the local extra faces for the interpolation. Note: could
        // let demand-driven lduAddr() trigger it but just to make sure.
        updateAddressing();

        // Addressing and/or weights have changed. Make interpolated cells
        // up to date with donors
        interpolateFields();

        return true;
    }

    return false;
}


Foam::word Foam::dynamicOversetFvMesh::baseName(const word& name)
{
    if (name.endsWith("_0"))
    {
        return baseName(name.substr(0, name.size()-2));
    }
    else
    {
        return name;
    }
}


bool Foam::dynamicOversetFvMesh::interpolateFields()
{
    // Add the stencil suppression list
    wordHashSet suppressed(Stencil::New(*this).nonInterpolatedFields());

    // Use whatever the solver has set up as suppression list
    const dictionary* dictPtr
    (
        this->schemesDict().findDict("oversetInterpolationSuppressed")
    );
    if (dictPtr)
    {
        suppressed.insert(dictPtr->toc());
    }

    interpolate<volScalarField>(suppressed);
    interpolate<volVectorField>(suppressed);
    interpolate<volSphericalTensorField>(suppressed);
    interpolate<volSymmTensorField>(suppressed);
    interpolate<volTensorField>(suppressed);

    return true;
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
            dimensionedScalar(dimless, Zero),
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
            dimensionedScalar(dimless, Zero),
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
    if (debug)
    {
        const cellCellStencilObject& overlap = Stencil::New(*this);
        const labelIOList& zoneID = overlap.zoneID();
        const labelListList& cellStencil = overlap.cellStencil();

        labelList donorZoneID(zoneID);
        overlap.cellInterpolationMap().distribute(donorZoneID);

        forAll(cellStencil, cellI)
        {
            const labelList& stencil = cellStencil[cellI];
            if (stencil.size())
            {
                donorZoneID[cellI] = zoneID[stencil[0]];
                for (label i = 1; i < stencil.size(); i++)
                {
                    if (zoneID[stencil[i]] != donorZoneID[cellI])
                    {
                        WarningInFunction << "Mixed donor meshes for cell "
                            << cellI << " at " << C()[cellI]
                            << " donors:" << UIndirectList<point>(C(), stencil)
                            << endl;
                        donorZoneID[cellI] = -2;
                    }
                }
            }
        }

        volScalarField volDonorZoneID
        (
            IOobject
            (
                "donorZoneID",
                this->time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            *this,
            dimensionedScalar("minOne", dimless, scalar(-1)),
            zeroGradientFvPatchScalarField::typeName
        );
        forAll(donorZoneID, celli)
        {
            volDonorZoneID[celli] = donorZoneID[celli];
        }
        //- Do not correctBoundaryConditions since re-interpolates!
        //volDonorZoneID.correctBoundaryConditions();
        volDonorZoneID.writeObject(fmt, ver, cmp, valid);
    }

    return ok;
}


// ************************************************************************* //
