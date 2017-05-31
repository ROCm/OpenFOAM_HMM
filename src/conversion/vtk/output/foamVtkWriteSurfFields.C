/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "foamVtkWriteSurfFields.H"
#include "OFstream.H"
#include "emptyFvsPatchFields.H"
#include "fvsPatchFields.H"
#include "surfaceFields.H"
#include "foamVtkOutput.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::foamVtkOutput::writeSurfFields
(
    const fvMesh& mesh,
    const fileName& baseName,
    const foamVtkOutput::outputOptions outOpts,
    const UPtrList<const surfaceVectorField>& surfVectorFields
)
{
    outputOptions opts(outOpts);
    opts.legacy(true);  // Legacy only, no append

    std::ofstream os((baseName + (opts.legacy() ? ".vtk" : ".vtp")).c_str());
    autoPtr<foamVtkOutput::formatter> format = opts.newFormatter(os);

    if (opts.legacy())
    {
        foamVtkOutput::legacy::fileHeader(format(), "surfaceFields")
            << "DATASET POLYDATA" << nl;
    }

    const pointField& fc = mesh.faceCentres();

    os << "POINTS " << mesh.nFaces() << " float" << nl;

    foamVtkOutput::writeList(format(), fc);
    format().flush();

    foamVtkOutput::legacy::pointDataHeader
    (
        os,
        mesh.nFaces(),
        surfVectorFields.size()
    );


    // surfVectorFields
    forAll(surfVectorFields, fieldi)
    {
        const surfaceVectorField& svf = surfVectorFields[fieldi];

        os  << svf.name() << " 3 " << mesh.nFaces() << " float" << nl;
        for (label facei=0; facei < mesh.nInternalFaces(); ++facei)
        {
            foamVtkOutput::write(format(), svf[facei]);
        }

        forAll(svf.boundaryField(), patchi)
        {
            const fvPatch& pp = mesh.boundary()[patchi];
            const fvsPatchVectorField& pf = svf.boundaryField()[patchi];

            if (isA<emptyFvsPatchVectorField>(pf))
            {
                // Note: loop over polypatch size, not fvpatch size.
                forAll(pp.patch(), i)
                {
                    foamVtkOutput::write(format(), vector::zero);
                }
            }
            else
            {
                foamVtkOutput::writeList(format(), pf);
            }
        }

        format().flush();
    }
}


// ************************************************************************* //
