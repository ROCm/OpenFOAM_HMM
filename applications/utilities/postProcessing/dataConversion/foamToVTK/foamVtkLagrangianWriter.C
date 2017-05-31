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
    const word& cloudName,
    const fileName& baseName,
    const foamVtkOutput::outputOptions outOpts,
    const bool dummyCloud
)
:
    mesh_(mesh),
    format_(),
    cloudName_(cloudName),
    os_(),
    nParcels_(0)
{
    outputOptions opts(outOpts);
    opts.legacy(true);  // Legacy only, no append

    os_.open((baseName + (opts.legacy() ? ".vtk" : ".vtp")).c_str());

    format_ = opts.newFormatter(os_);

    if (opts.legacy())
    {
        foamVtkOutput::legacy::fileHeader(format(), mesh_.time().caseName())
            << "DATASET POLYDATA" << nl;
    }

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


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::foamVtkOutput::lagrangianWriter::~lagrangianWriter()
{}


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
