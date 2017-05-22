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

#include "foamVtkLagrangianWriter.H"
#include "Cloud.H"
#include "passiveParticle.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::foamVtkOutput::lagrangianWriter::lagrangianWriter
(
    const fvMesh& mesh,
    const bool binary,
    const fileName& fName,
    const word& cloudName,
    const bool dummyCloud
)
:
    mesh_(mesh),
    format_(),
    cloudName_(cloudName),
    os_(fName.c_str()),
    nParcels_(0)
{
    format_ = foamVtkOutput::newFormatter
    (
        os_,
        (
            binary
          ? foamVtkOutput::LEGACY_BINARY
          : foamVtkOutput::LEGACY_ASCII
        )
    );

    foamVtkOutput::legacy::fileHeader(format(), mesh_.time().caseName())
        << "DATASET POLYDATA" << nl;

    if (dummyCloud)
    {
        os_ << "POINTS " << nParcels_ << " float" << nl;
    }
    else
    {
        Cloud<passiveParticle> parcels(mesh, cloudName_, false);

        nParcels_ = parcels.size();

        os_ << "POINTS " << nParcels_ << " float" << nl;

        forAllConstIters(parcels, iter)
        {
            const point& pt = iter().position();

            foamVtkOutput::write(format(), pt);
        }
        format().flush();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::foamVtkOutput::lagrangianWriter::beginParcelData
(
    const label nFields
)
{
    foamVtkOutput::legacy::pointDataHeader(os_, nParcels_, nFields);
}


void Foam::foamVtkOutput::lagrangianWriter::endParcelData()
{}


// ************************************************************************* //
