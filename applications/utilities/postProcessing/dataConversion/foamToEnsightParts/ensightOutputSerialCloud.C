/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "ensightOutputSerialCloud.H"
#include "ensightPTraits.H"

#include "passiveParticle.H"
#include "IOField.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::ensightSerialCloud::writePositions
(
    const polyMesh& mesh,
    const word& cloudName,
    autoPtr<ensightFile> output
)
{
    label nTotParcels = 0;
    autoPtr<Cloud<passiveParticle>> cloudPtr;

    cloudPtr.reset(new Cloud<passiveParticle>(mesh, cloudName, false));
    nTotParcels = cloudPtr().size();

    Cloud<passiveParticle> parcels(mesh, cloudName, false);

    if (Pstream::master())
    {
        ensightFile& os = output();
        os.beginParticleCoordinates(nTotParcels);

        // binary write is Ensight6 - first ids, then positions
        if (os.format() == IOstream::BINARY)
        {
            // 1-index
            for (label parcelId = 0; parcelId < nTotParcels; ++parcelId)
            {
                os.write(parcelId+1);
            }


            forAllConstIter(Cloud<passiveParticle>, cloudPtr(), elmnt)
            {
                const vector p(elmnt().position());

                os.write(p.x());
                os.write(p.y());
                os.write(p.z());
            }
        }
        else
        {
            // ASCII id + position together

            label parcelId = 0;
            forAllConstIter(Cloud<passiveParticle>, cloudPtr(), elmnt)
            {
                const vector p(elmnt().position());

                os.write(++parcelId, 8); // unusual width
                os.write(p.x());
                os.write(p.y());
                os.write(p.z());
                os.newline();
            }
        }
    }
}


// ************************************************************************* //
