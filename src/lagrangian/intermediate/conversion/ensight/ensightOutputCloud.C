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

\*---------------------------------------------------------------------------*/

#include "ensightOutputCloud.H"
#include "fvMesh.H"
#include "Cloud.H"
#include "passiveParticle.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    //- Binary output
    static inline void writeMeasured
    (
        ensightFile& os,
        const pointField& points
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
    static inline label writeMeasured
    (
        ensightFile& os,
        label pointId,
        const pointField& points
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
    autoPtr<ensightFile>& output,
    Pstream::commsTypes comm
)
{
    pointField positions;

    if (exists)
    {
        Cloud<passiveParticle> parcels(mesh, cloudName, false);

        positions.resize(parcels.size());

        auto outIter = positions.begin();

        for (const passiveParticle& p : parcels)
        {
            *outIter = p.position();
            ++outIter;
        }
    }


    // Total number of parcels on all processes
    const label nTotParcels = returnReduce(positions.size(), sumOp<label>());

    // Update the exists/not exists information (for return value)
    exists = nTotParcels;

    if (Pstream::master())
    {
        ensightFile& os = output();

        os.beginParticleCoordinates(nTotParcels);
        if (!exists)
        {
            return exists;  // DONE
        }

        if (os.format() == IOstream::BINARY)
        {
            // binary write is Ensight6 - first ids, then positions

            // 1-index
            for (label parcelId = 1; parcelId <= nTotParcels; ++parcelId)
            {
                os.write(parcelId);
            }

            // Master
            writeMeasured(os, positions);

            // Slaves
            for (const int slave : Pstream::subProcs())
            {
                IPstream fromSlave(comm, slave);
                pointField recv(fromSlave);

                writeMeasured(os, recv);
            }
        }
        else
        {
            // ASCII id + position together
            label parcelId = 0;

            // Master
            parcelId = writeMeasured(os, parcelId, positions);

            // Slaves
            for (const int slave : Pstream::subProcs())
            {
                IPstream fromSlave(comm, slave);
                pointField recv(fromSlave);

                parcelId = writeMeasured(os, parcelId, recv);
            }
        }
    }
    else if (nTotParcels)
    {
        // SLAVE, and data exist
        OPstream toMaster(comm, Pstream::masterNo());

        toMaster
            << positions;
    }

    return exists;
}


// ************************************************************************* //
