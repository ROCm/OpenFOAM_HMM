/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 DLR
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "zoneDistribute.H"
#include "DynamicField.H"
#include "syncTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<typename Type>
Type Foam::zoneDistribute::getLocalValue
(
    const GeometricField<Type, fvPatchField, volMesh>& phi,
    const label localIdx
) const
{
    if (localIdx < mesh_.nCells()) // internal: cellI
    {
        return phi[localIdx];
    }

    return faceValue(phi,localIdx);
}


template<typename Type>
Type Foam::zoneDistribute::faceValue
(
    const GeometricField<Type, fvPatchField, volMesh>& phi,
    const label localIdx
) const
{
    const label faceI = localIdx + mesh_.nInternalFaces() - mesh_.nCells();

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    // Boundary face. Find out which face of which patch
    const label patchI = pbm.whichPatch(faceI);

    if (patchI < 0 || patchI >= pbm.size())
    {
       FatalErrorInFunction
            << "Cannot find patch for face " << faceI
            << abort(FatalError);
    }

    const polyPatch& pp = pbm[patchI];

    const label patchFaceI = pp.whichFace(faceI);

    return phi.boundaryField()[patchI][patchFaceI];
}


template<typename Type>
Type Foam::zoneDistribute::getValue
(
    const GeometricField<Type, fvPatchField, volMesh>& phi,
    const Map<Type>& valuesFromOtherProc,
    const label gblIdx
) const
{
    if (globalNumbering_.isLocal(gblIdx))
    {
        const label localIdx = globalNumbering_.toLocal(gblIdx);
        return getLocalValue(phi,localIdx);
    }
    else // from other proc
    {
        return valuesFromOtherProc[gblIdx];
    }
}


template<typename Type>
Foam::Map<Foam::Field<Type>> Foam::zoneDistribute::getFields
(
    const boolList& zone,
    const GeometricField<Type, fvPatchField, volMesh>& phi
)
{
    if (zone.size() != phi.size())
    {
        FatalErrorInFunction
            << "size of zone: " << zone.size()
            << "size of phi:" <<  phi.size()
            << "do not match. Did the mesh change?"
            << exit(FatalError);
    }


    // Get values from other proc
    Map<Type> neiValues = getDatafromOtherProc(zone, phi);

    Map<Field<Type>> stencilWithValues;

    DynamicField<Type> tmpField(128);

    forAll(zone, celli)
    {
        if (zone[celli])
        {
            tmpField.clear();

            for (const label gblIdx : stencil_[celli])
            {
                tmpField.append(getValue(phi,neiValues,gblIdx));
            }

            stencilWithValues.emplace(celli, tmpField);
        }
    }

    return stencilWithValues;
}


template<typename Type>
Foam::Map<Type> Foam::zoneDistribute::getDatafromOtherProc
(
    const boolList& zone,
    const GeometricField<Type, fvPatchField, volMesh>& phi
)
{
    if (zone.size() != phi.size())
    {
        FatalErrorInFunction
            << "size of zone: " << zone.size()
            << "size of phi:" <<  phi.size()
            << "do not match. Did the mesh change?"
            << exit(FatalError);
    }


    // Get values from other proc
    Map<Type> neiValues;

    if (UPstream::parRun())
    {
        const bool commsPrepared = (sendTo_.size() == UPstream::nProcs());

        if (!commsPrepared)
        {
            FatalErrorInFunction
                << "The send/recv information not initialized - "
                << "likely that setUpCommforZone() was not called"
                << endl;

            // But don't exit/abort for now
        }


        List<Map<Type>> sendValues(UPstream::nProcs());

        forAll(send_, proci)
        {
            for (const label sendIdx : send_[proci])
            {
                sendValues[proci].insert
                (
                    sendIdx,
                    getLocalValue(phi, globalNumbering_.toLocal(sendIdx))
                );
            }
        }


        PstreamBuffers pBufs(UPstream::commsTypes::nonBlocking);
        labelList recvSizes;

        // Stream data into buffer
        for (const int proci : UPstream::allProcs())
        {
            if (proci != UPstream::myProcNo() && !sendValues[proci].empty())
            {
                // Put data into send buffer
                UOPstream toProc(proci, pBufs);
                toProc << sendValues[proci];
            }
        }


        if (commsPrepared)
        {
            DynamicList<label> sendProcs;
            DynamicList<label> recvProcs;

            for (const int proci : UPstream::allProcs())
            {
                if (sendTo_[proci])
                {
                    sendProcs.append(proci);
                }
                if (recvFrom_[proci])
                {
                    recvProcs.append(proci);
                }
            }

            pBufs.finishedSends(sendProcs, recvProcs, recvSizes);
        }
        else
        {
            // Fallback to all-to-all communication

            pBufs.finishedSends(recvSizes);
        }


        for (const int proci : UPstream::allProcs())
        {
            // Do not use recvFrom_[proci] here
            // - could be incorrect if comms are not prepared!

            if
            (
                proci != UPstream::myProcNo()
             && (recvSizes[proci] > 0)
            )
            {
                UIPstream fromProc(proci, pBufs);
                Map<Type> tmpValue(fromProc);

                neiValues += tmpValue;
            }
        }
    }

    return neiValues;
}


// ************************************************************************* //
