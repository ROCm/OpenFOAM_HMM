/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "ensightOutputCloud.H"
#include "fvMesh.H"
#include "Cloud.H"
#include "passiveParticle.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    //- Binary output
    static inline void writeMeasured_binary
    (
        ensightFile& os,
        const UList<point>& points
    )
    {
        for (const point& p : points)
        {
            os.write(p.x());
            os.write(p.y());
            os.write(p.z());
        }
    }

    //- ASCII output. Id + position together
    static inline label writeMeasured_ascii
    (
        ensightFile& os,
        label pointId,
        const UList<point>& points
    )
    {
        for (const point& p : points)
        {
            os.write(++pointId, 8); // 1-index and an unusual width
            os.write(p.x());
            os.write(p.y());
            os.write(p.z());
            os.newline();
        }

        return pointId;
    }

} // End namespace Foam


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

bool Foam::ensightOutput::writeCloudPositions
(
    const fvMesh& mesh,
    const word& cloudName,
    bool exists,
    autoPtr<ensightFile>& output
)
{
    label nLocalParcels(0);
    autoPtr<Cloud<passiveParticle>> parcelsPtr;

    if (exists)
    {
        parcelsPtr.reset(new Cloud<passiveParticle>(mesh, cloudName, false));
        nLocalParcels = parcelsPtr().size();
    }

    // Total number of parcels on all processes
    const label nTotParcels = returnReduce(nLocalParcels, sumOp<label>());

    if (Pstream::master())
    {
        ensightFile& os = output();
        os.beginParticleCoordinates(nTotParcels);
    }

    if (!nTotParcels)
    {
        return false;  // DONE
    }


    // Size information (offsets are irrelevant)
    const globalIndex procAddr
    (
        UPstream::listGatherValues<label>(nLocalParcels),
        globalIndex::SIZES
    );


    DynamicList<point> positions;
    positions.reserve(Pstream::master() ? procAddr.maxSize() : nLocalParcels);

    // Extract positions
    if (parcelsPtr)
    {
        const auto& parcels = *parcelsPtr;

        positions.resize_nocopy(parcels.size());  // same as nLocalParcels

        auto outIter = positions.begin();

        for (const passiveParticle& p : parcels)
        {
            *outIter = p.position();
            ++outIter;
        }

        parcelsPtr.reset(nullptr);
    }

    if (Pstream::master())
    {
        ensightFile& os = output();
        const bool isBinaryOutput = (os.format() == IOstream::BINARY);

        label parcelId = 0;

        if (isBinaryOutput)
        {
            // NB: binary write is Ensight6 - first ids, then positions

            // 1-index
            for (label id = 1; id <= nTotParcels; ++id)
            {
                os.write(id);
            }

            // Write master data
            writeMeasured_binary(os, positions);
        }
        else
        {
            // NB: ascii write is (id + position) together

            // Write master data
            parcelId = writeMeasured_ascii(os, parcelId, positions);
        }


        // Receive and write
        for (const label proci : procAddr.subProcs())
        {
            positions.resize_nocopy(procAddr.localSize(proci));
            UIPstream::read
            (
                UPstream::commsTypes::scheduled,
                proci,
                positions.data_bytes(),
                positions.size_bytes()
            );

            if (isBinaryOutput)
            {
                writeMeasured_binary(os, positions);
            }
            else
            {
                parcelId = writeMeasured_ascii(os, parcelId, positions);
            }
        }
    }
    else
    {
        // Send
        UOPstream::write
        (
            UPstream::commsTypes::scheduled,
            Pstream::masterNo(),
            positions.cdata_bytes(),
            positions.size_bytes()
        );
    }

    return true;
}


// ************************************************************************* //
