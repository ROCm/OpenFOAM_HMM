/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "vtkWrite.H"
#include "dictionary.H"
#include "Time.H"
#include "areaFields.H"
#include "foamVtkInternalWriter.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(vtkWrite, 0);
    addToRunTimeSelectionTable(functionObject, vtkWrite, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::vtkWrite::vtkWrite
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeOpts_(vtk::formatType::INLINE_BASE64),
    selectFields_(),
    dirName_("VTK")
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::vtkWrite::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    //
    // writer options - default is xml base64
    //
    writeOpts_ = vtk::formatType::INLINE_BASE64;
    if (dict.lookupOrDefault<bool>("legacy", false))
    {
        writeOpts_.legacy(true);
    }

    writeOpts_.ascii
    (
        dict.found("format")
     && (IOstream::formatEnum(dict.lookup("format")) == IOstream::ASCII)
    );

    // FUTURE?
    // writeOpts_.precision
    // (
    //     dict.lookupOrDefault
    //     (
    //         "writePrecision",
    //         IOstream::defaultPrecision()
    //     )
    // );

    // Info<< type() << " " << name() << " output-format: "
    //     << writeOpts_.description() << nl;

    //
    // other options
    //
    dict.readIfPresent("directory", dirName_);

    writeIds_ = dict.lookupOrDefault<bool>("writeIds", false);


    //
    // output fields
    //
    dict.lookup("fields") >> selectFields_;
    wordRes::inplaceUniq(selectFields_);

    return true;
}


bool Foam::functionObjects::vtkWrite::execute()
{
    return true;
}


bool Foam::functionObjects::vtkWrite::write()
{
    // const word timeDesc =
    //     useTimeName ? time_.timeName() : Foam::name(time_.timeIndex());

    const word timeDesc = time_.timeName();

    fileName vtkDir = dirName_;
    if (!vtkDir.isAbsolute())
    {
        vtkDir = time_.path()/vtkDir;
    }
    mkDir(vtkDir);

    string vtkName = time_.caseName();

    if (Pstream::parRun())
    {
        // Strip off leading casename, leaving just processor_DDD ending.
        string::size_type i = vtkName.rfind("processor");

        if (i != string::npos)
        {
            vtkName = vtkName.substr(i);
        }
    }

    // internal mesh
    {
        const fileName outputName
        (
            vtkDir/vtkName
          + "_"
          + timeDesc
        );

        Info<< name() << " output Time: " << time_.timeName() << nl
            << "    Internal  : " << outputName << endl;

        // Number of fields to be written: only needed for legacy vtk format
        label nVolFields = 0;
        if (writeOpts_.legacy())
        {
            nVolFields =
            (
                (writeIds_ ? 1 : 0)
              + countFields<volScalarField>()
              + countFields<volVectorField>()
              + countFields<volSphericalTensorField>()
              + countFields<volSymmTensorField>()
              + countFields<volTensorField>()
            );
        }

        vtk::vtuCells vtuMeshCells
        (
            mesh_,
            writeOpts_,
            true  // decompose
        );

        // Write mesh
        vtk::internalWriter writer
        (
            mesh_,
            vtuMeshCells,
            outputName,
            writeOpts_
        );

        // CellData
        {
            writer.beginCellData(nVolFields);

            // Write cellID field
            if (writeIds_)
            {
                writer.writeCellIDs();
            }

            // Write volFields
            writeFields<volScalarField>(writer);
            writeFields<volVectorField>(writer);
            writeFields<volSphericalTensorField>(writer);
            writeFields<volSymmTensorField>(writer);
            writeFields<volTensorField>(writer);

            writer.endCellData();
        }

        writer.writeFooter();
    }

    return true;
}


// ************************************************************************* //
