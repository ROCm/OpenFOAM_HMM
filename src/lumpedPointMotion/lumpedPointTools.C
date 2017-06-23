/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "lumpedPointTools.H"
#include "IFstream.H"
#include "IOobjectList.H"
#include "volFields.H"

#include "lumpedPointDisplacementPointPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // file-scope
    template<class GeoFieldType>
    static autoPtr<GeoFieldType> loadPointField
    (
        const pointMesh::Mesh& mesh,
        const IOobject* io
    )
    {
        if (io && io->headerClassName() == GeoFieldType::typeName)
        {
            Info<< "Reading " << GeoFieldType::typeName
                << ' ' << io->name() << endl;

            return autoPtr<GeoFieldType>
            (
                new GeoFieldType
                (
                    IOobject
                    (
                        io->name(),
                        io->instance(),
                        io->local(),
                        io->db(),
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE,
                        io->registerObject()
                    ),
                    mesh
                )
            );
        }

        return autoPtr<GeoFieldType>();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::List<Foam::lumpedPointStateTuple>
Foam::lumpedPointTools::lumpedPointStates(Istream& is)
{
    dictionary contents(is);
    List<dictionary> entries(contents.lookup("response"));

    DynamicList<Tuple2<scalar, lumpedPointState>> states(entries.size());

    for (const dictionary& dict : entries)
    {
        states.append
        (
            lumpedPointStateTuple
            (
                readScalar(dict.lookup("time")),
                lumpedPointState(dict)
            )
        );
    }

    return states.shrink();
}


Foam::List<Foam::lumpedPointStateTuple>
Foam::lumpedPointTools::lumpedPointStates(const fileName& file)
{
    IFstream is(file);
    return lumpedPointStates(is);
}


Foam::pointIOField
Foam::lumpedPointTools::points0Field(const polyMesh& mesh)
{
    pointIOField pts
    (
        IOobject
        (
            "points",
            mesh.time().constant(),
            polyMesh::meshSubDir,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false // Do not re-register
        )
    );

    return pts;
}



Foam::labelList
Foam::lumpedPointTools::lumpedPointPatchList(const pointVectorField& pvf)
{
    return lumpedPointDisplacementPointPatchVectorField::patchIds(pvf);
}


Foam::labelList Foam::lumpedPointTools::lumpedPointPatchList
(
    const polyMesh& mesh
)
{
    IOobjectList objects0(mesh, "0");

    pointMesh pMesh(mesh);

    autoPtr<pointVectorField> displacePtr = loadPointField<pointVectorField>
    (
        pMesh,
        objects0.lookup("pointDisplacement")
    );

    if (!displacePtr.valid())
    {
        Info<< "no valid pointDisplacement" << endl;
        return labelList();
    }

    return lumpedPointPatchList(displacePtr());
}


// ************************************************************************* //
