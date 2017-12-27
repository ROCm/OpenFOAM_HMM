/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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
#include "passiveParticle.H"
#include "Cloud.H"
#include "pointList.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::ensightCloud::writePositions
(
    const fvMesh& mesh,
    const word& cloudName,
    const bool exists,
    autoPtr<ensightFile>& output
)
{
    // Total number of parcels on all processes
    label nTotParcels = 0;
    autoPtr<Cloud<passiveParticle>> cloudPtr;

    if (exists)
    {
        cloudPtr.reset(new Cloud<passiveParticle>(mesh, cloudName, false));
        nTotParcels = cloudPtr().size();
    }
    reduce(nTotParcels, sumOp<label>());

    if (Pstream::master())
    {
        ensightFile& os = output();

        os.beginParticleCoordinates(nTotParcels);
        if (!nTotParcels)
        {
            return;  // DONE
        }

        if (os.format() == IOstream::BINARY)
        {
            // binary write is Ensight6 - first ids, then positions

            // 1-index
            for (label parcelId = 0; parcelId < nTotParcels; ++parcelId)
            {
                os.write(parcelId+1);
            }

            // Master
            forAllConstIter(Cloud<passiveParticle>, cloudPtr(), elmnt)
            {
                const point p(elmnt().position());

                os.write(p.x());
                os.write(p.y());
                os.write(p.z());
            }

            // Slaves
            for (int slave=1; slave<Pstream::nProcs(); ++slave)
            {
                IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
                pointList points(fromSlave);

                forAll(points, pti)
                {
                    const point& p = points[pti];

                    os.write(p.x());
                    os.write(p.y());
                    os.write(p.z());
                }
            }
        }
        else
        {
            // ASCII id + position together

            label parcelId = 0;
            forAllConstIter(Cloud<passiveParticle>, cloudPtr(), elmnt)
            {
                const point p(elmnt().position());

                os.write(++parcelId, 8); // unusual width
                os.write(p.x());
                os.write(p.y());
                os.write(p.z());
                os.newline();
            }

            // Slaves
            for (int slave=1; slave<Pstream::nProcs(); ++slave)
            {
                IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
                pointList points(fromSlave);

                forAll(points, pti)
                {
                    const point& p = points[pti];

                    os.write(++parcelId, 8); // unusual width
                    os.write(p.x());
                    os.write(p.y());
                    os.write(p.z());
                    os.newline();
                }
            }
        }
    }
    else if (nTotParcels)
    {
        // SLAVE, and data exist
        pointList points(cloudPtr().size());

        label pti = 0;
        forAllConstIter(Cloud<passiveParticle>, cloudPtr(), elmnt)
        {
            const point p(elmnt().position());
            points[pti++] = p;
        }

        {
            OPstream toMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );

            toMaster
                << points;
        }
    }
}


// ************************************************************************* //
