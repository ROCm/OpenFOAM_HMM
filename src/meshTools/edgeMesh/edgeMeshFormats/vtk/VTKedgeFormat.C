/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "VTKedgeFormat.H"
#include "Fstream.H"
#include "clock.H"
#include "vtkUnstructuredReader.H"
#include "Time.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fileFormats::VTKedgeFormat::writeHeader
(
    Ostream& os,
    const pointField& pointLst
)
{
    // Write header
    os  << "# vtk DataFile Version 2.0" << nl
        << "featureEdgeMesh written " << clock::dateTime().c_str() << nl
        << "ASCII" << nl
        << nl
        << "DATASET POLYDATA" << nl;

    // Write vertex coords
    os  << "POINTS " << pointLst.size() << " double" << nl;
    for (const point& pt : pointLst)
    {
        os  << float(pt.x()) << ' '
            << float(pt.y()) << ' '
            << float(pt.z()) << nl;
    }
}


void Foam::fileFormats::VTKedgeFormat::writeEdges
(
    Ostream& os,
    const UList<edge>& edgeLst
)
{
    os  << "LINES " << edgeLst.size() << ' ' << 3*edgeLst.size() << nl;

    for (const edge& e : edgeLst)
    {
        os  << "2 " << e[0] << ' ' << e[1] << nl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::VTKedgeFormat::VTKedgeFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileFormats::VTKedgeFormat::read
(
    const fileName& filename
)
{
    IFstream is(filename);
    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << filename
            << exit(FatalError);
    }

    // Use dummy Time for objectRegistry
    autoPtr<Time> dummyTimePtr(Time::New());

    objectRegistry obr
    (
        IOobject
        (
            "vtk::edgeFormat",
            *dummyTimePtr,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    // Construct reader to read file
    vtkUnstructuredReader reader(obr, is);


    // Extract lines
    storedPoints().transfer(reader.points());

    label nEdges = 0;
    forAll(reader.lines(), lineI)
    {
        nEdges += reader.lines()[lineI].size()-1;
    }
    storedEdges().setSize(nEdges);

    nEdges = 0;
    forAll(reader.lines(), lineI)
    {
        const labelList& verts = reader.lines()[lineI];
        for (label i = 1; i < verts.size(); i++)
        {
            storedEdges()[nEdges++] = edge(verts[i-1], verts[i]);
        }
    }

    return true;
}


void Foam::fileFormats::VTKedgeFormat::write
(
    const fileName& filename,
    const edgeMesh& eMesh
)
{
    OFstream os(filename);
    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot open file for writing " << filename
            << exit(FatalError);
    }

    writeHeader(os, eMesh.points());
    writeEdges(os, eMesh.edges());
}


// ************************************************************************* //
