/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 DHI
    Copyright (C) 2017-2020 OpenCFD Ltd.
    Copyright (c) 2017-2020, German Aerospace Center (DLR)
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

Application
    setAlphaField

Description
    Uses cutCellIso to create a volume fraction field from either a cylinder,
    a sphere or a plane.

    Original code supplied by Johan Roenby, DHI (2016)
    Modification Henning Scheufler, DLR (2019)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "triSurface.H"
#include "triSurfaceTools.H"

#include "implicitFunction.H"

#include "cutCellIso.H"
#include "OBJstream.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void isoFacesToFile
(
    const DynamicList<List<point>>& faces,
    const word& fileName
)
{
    // Writing isofaces to OBJ file for inspection in paraview
    OBJstream os(fileName + ".obj");

    if (Pstream::parRun())
    {
        // Collect points from all the processors
        List<DynamicList<List<point>>> allProcFaces(Pstream::nProcs());
        allProcFaces[Pstream::myProcNo()] = faces;
        Pstream::gatherList(allProcFaces);

        if (Pstream::master())
        {
            Info<< "Writing file: " << fileName << endl;

            for (const DynamicList<List<point>>& procFaces : allProcFaces)
            {
                for (const List<point>& facePts : procFaces)
                {
                    os.write(face(identity(facePts.size())), facePts);
                }
            }
        }
    }
    else
    {
        Info<< "Writing file: " << fileName << endl;

        for (const List<point>& facePts : faces)
        {
            os.write(face(identity(facePts.size())), facePts);
        }
    }
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Uses cutCellIso to create a volume fraction field from an "
        "implicit function."
    );

    #include "addDictOption.H"
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    const word dictName("setAlphaFieldDict");
    #include "setSystemMeshDictionaryIO.H"

    IOdictionary setAlphaFieldDict(dictIO);

    Info<< "Reading " << setAlphaFieldDict.name() << endl;

    const word fieldName = setAlphaFieldDict.get<word>("field");
    const bool invert = setAlphaFieldDict.getOrDefault("invert", false);
    const bool writeOBJ = setAlphaFieldDict.getOrDefault("writeOBJ", false);

    Info<< "Reading field " << fieldName << nl << endl;
    volScalarField alpha1
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    autoPtr<implicitFunction> func = implicitFunction::New
    (
        setAlphaFieldDict.get<word>("type"),
        setAlphaFieldDict
    );

    scalarField f(mesh.nPoints(), Zero);

    forAll(f, pi)
    {
        f[pi] = func->value(mesh.points()[pi]);
    };

    cutCellIso cutCell(mesh, f);

    DynamicList<List<point>> facePts;

    DynamicList<triSurface> surface;

    surfaceScalarField cellToCellDist(mag(mesh.delta()));

    forAll(alpha1, cellI)
    {
        label cellStatus = cutCell.calcSubCell(cellI, 0.0);

        if (cellStatus == -1)
        {
            alpha1[cellI] = 1;
        }
        else if (cellStatus == 1)
        {
            alpha1[cellI] = 0;
        }
        else if (cellStatus == 0)
        {
            if (mag(cutCell.faceArea()) != 0)
            {
                alpha1[cellI] = max(min(cutCell.VolumeOfFluid(), 1), 0);
                if (writeOBJ && (mag(cutCell.faceArea()) >= 1e-14))
                {
                    facePts.append(cutCell.facePoints());
                }
            }
        }
    }

    if (writeOBJ)
    {
        isoFacesToFile(facePts, fieldName + "0");
    }

    ISstream::defaultPrecision(18);

    if (invert)
    {
        alpha1 = scalar(1) - alpha1;
    }

    alpha1.correctBoundaryConditions();

    Info<< "Writing new alpha field " << alpha1.name() << endl;
    alpha1.write();

    const scalarField& alpha = alpha1.internalField();

    Info<< "sum(alpha*V):" << gSum(mesh.V()*alpha)
        << ", 1-max(alpha1): " << 1 - gMax(alpha)
        << " min(alpha1): " << gMin(alpha) << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
