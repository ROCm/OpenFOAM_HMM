/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

#include "vtkCloud.H"
#include "Cloud.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "foamVtkOutputOptions.H"
#include "addToRunTimeSelectionTable.H"
#include "pointList.H"
#include "stringOps.H"
#include <fstream>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(vtkCloud, 0);
    addToRunTimeSelectionTable(functionObject, vtkCloud, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::vtkCloud::writeVerts
(
    autoPtr<vtk::formatter>& format,
    const label nParcels
) const
{
    if (Pstream::master())
    {
        format().tag(vtk::fileTag::VERTS);

        // Same payload throughout
        const uint64_t payLoad = (nParcels * sizeof(label));

        //
        // 'connectivity'
        // = linear mapping onto points
        //
        {
            format().openDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY)
                .closeTag();

            format().writeSize(payLoad);
            for (label i=0; i < nParcels; ++i)
            {
                format().write(i);
            }
            format().flush();

            format().endDataArray();
        }

        //
        // 'offsets'  (connectivity offsets)
        // = linear mapping onto points (with 1 offset)
        //
        {
            format().openDataArray<label>(vtk::dataArrayAttr::OFFSETS)
                .closeTag();

            format().writeSize(payLoad);
            for (label i=0; i < nParcels; ++i)
            {
                format().write(i+1);
            }
            format().flush();

            format().endDataArray();
        }

        format().endTag(vtk::fileTag::VERTS);
    }
}


bool Foam::functionObjects::vtkCloud::writeCloud
(
    const fileName& outputName,
    const word& cloudName
)
{
    const auto* objPtr = mesh_.findObject<cloud>(cloudName);
    if (!objPtr)
    {
        return false;
    }

    objectRegistry obrTmp
    (
        IOobject
        (
            "vtk::vtkCloud::" + cloudName,
            mesh_.time().constant(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    objPtr->writeObjects(obrTmp);

    const auto* pointsPtr = obrTmp.findObject<vectorField>("position");

    if (!pointsPtr)
    {
        // This should be impossible
        return false;
    }

    // Total number of parcels on all processes
    label nTotParcels = pointsPtr->size();
    reduce(nTotParcels, sumOp<label>());

    if (pruneEmpty_ && !nTotParcels)
    {
        return false;
    }

    std::ofstream os;
    autoPtr<vtk::formatter> format;

    // Header
    if (Pstream::master())
    {
        os.open(outputName);
        format = writeOpts_.newFormatter(os);

        // XML (inline)
        format()
            .xmlHeader()
            .xmlComment
            (
                "cloud=" + cloudName
              + " time=" + time_.timeName()
              + " index=" + Foam::name(time_.timeIndex())
            )
            .beginVTKFile(vtk::fileTag::POLY_DATA, "0.1");

        // Begin piece
        if (useVerts_)
        {
            format()
                .openTag(vtk::fileTag::PIECE)
                .xmlAttr(vtk::fileAttr::NUMBER_OF_POINTS, nTotParcels)
                .xmlAttr(vtk::fileAttr::NUMBER_OF_VERTS, nTotParcels)
                .closeTag();
        }
        else
        {
            format()
                .openTag(vtk::fileTag::PIECE)
                .xmlAttr(vtk::fileAttr::NUMBER_OF_POINTS, nTotParcels)
                .closeTag();
        }
    }


    // Points
    if (Pstream::master())
    {
        const uint64_t payLoad = (nTotParcels * 3 * sizeof(float));

        format().tag(vtk::fileTag::POINTS)
            .openDataArray<float,3>(vtk::dataArrayAttr::POINTS)
            .closeTag();

        format().writeSize(payLoad);

        // Master
        vtk::writeList(format(), *pointsPtr);

        // Slaves - recv
        for (int slave=1; slave<Pstream::nProcs(); ++slave)
        {
            IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
            pointList points(fromSlave);

            vtk::writeList(format(), points);
        }

        format().flush();

        format()
            .endDataArray()
            .endTag(vtk::fileTag::POINTS);

        if (useVerts_)
        {
            writeVerts(format, nTotParcels);
        }
    }
    else
    {
        // Slaves - send

        OPstream toMaster
        (
            Pstream::commsTypes::scheduled,
            Pstream::masterNo()
        );

        toMaster
            << *pointsPtr;
    }


    // Prevent any possible conversion of positions as a field
    obrTmp.filterKeys
    (
        [](const word& k)
        {
            return k.startsWith("position")
                || k.startsWith("coordinate");
        },
        true  // prune
    );

    // Restrict to specified fields
    if (selectFields_.size())
    {
        obrTmp.filterKeys(selectFields_);
    }


    // Write fields
    const vtk::fileTag dataType =
    (
        useVerts_
      ? vtk::fileTag::CELL_DATA
      : vtk::fileTag::POINT_DATA
    );

    if (Pstream::master())
    {
        format().tag(dataType);
    }

    DynamicList<word> written(obrTmp.size());

    written.append(writeFields<label>(format, obrTmp, nTotParcels));
    written.append(writeFields<scalar>(format, obrTmp, nTotParcels));
    written.append(writeFields<vector>(format, obrTmp, nTotParcels));

    if (Pstream::master())
    {
        format().endTag(dataType);
    }

    // Footer
    if (Pstream::master())
    {
        // slight cheat. </Piece> too
        format().endTag(vtk::fileTag::PIECE);

        format().endTag(vtk::fileTag::POLY_DATA)
            .endVTKFile();
    }


    // Record information into the state (all processors)
    //
    // foName
    // {
    //     cloudName
    //     {
    //         file   "<case>/VTK/cloud1_000.vtp";
    //         fields (U T rho);
    //     }
    // }

    dictionary propsDict;

    // Use case-local filename and "<case>" shortcut for readable output
    // and possibly relocation of files

    fileName fName(outputName.relative(stringOps::expand("<case>")));
    if (fName.isAbsolute())
    {
        propsDict.add("file", fName);
    }
    else
    {
        propsDict.add("file", "<case>"/fName);
    }
    propsDict.add("fields", written);

    setObjectProperty(name(), cloudName, propsDict);

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::vtkCloud::vtkCloud
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeOpts_(vtk::formatType::INLINE_BASE64),
    printf_(),
    useVerts_(false),
    pruneEmpty_(false),
    selectClouds_(),
    selectFields_(),
    dirName_("VTK"),
    series_()
{
    if (postProcess)
    {
        // Disable for post-process mode.
        // Emit as FatalError for the try/catch in the caller.
        FatalError
            << type() << " disabled in post-process mode"
            << exit(FatalError);
    }

    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::vtkCloud::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    //
    // writer options - default is xml base64. Legacy is not desired.
    //
    writeOpts_ = vtk::formatType::INLINE_BASE64;

    writeOpts_.ascii
    (
        dict.found("format")
     && (IOstream::formatEnum(dict.get<word>("format")) == IOstream::ASCII)
    );

    writeOpts_.append(false);  // No append supported
    writeOpts_.legacy(false);  // No legacy supported

    writeOpts_.precision
    (
        dict.lookupOrDefault
        (
            "writePrecision",
            IOstream::defaultPrecision()
        )
    );

    // Info<< type() << " " << name() << " output-format: "
    //     << writeOpts_.description() << nl;

    int padWidth = dict.lookupOrDefault<int>("width", 8);

    // Appropriate printf format - Enforce min/max sanity limits
    if (padWidth < 1 || padWidth > 31)
    {
        printf_.clear();
    }
    else
    {
        printf_ = "%0" + std::to_string(padWidth) + "d";
    }

    // useTimeName_ = dict.lookupOrDefault<bool>("useTimeName", false);

    useVerts_ = dict.lookupOrDefault<bool>("cellData", false);
    pruneEmpty_ = dict.lookupOrDefault<bool>("prune", false);


    //
    // other options
    //
    dict.readIfPresent("directory", dirName_);

    selectClouds_.clear();
    dict.readIfPresent("clouds", selectClouds_);

    if (selectClouds_.empty())
    {
        selectClouds_.resize(1);
        selectClouds_.first() =
            dict.lookupOrDefault<word>("cloud", cloud::defaultName);
    }

    selectFields_.clear();
    dict.readIfPresent("fields", selectFields_);

    return true;
}


bool Foam::functionObjects::vtkCloud::execute()
{
    return true;
}


bool Foam::functionObjects::vtkCloud::write()
{
    const wordList cloudNames(mesh_.sortedNames<cloud>(selectClouds_));

    if (cloudNames.empty())
    {
        return true;  // skip - not available
    }

//     const word timeDesc =
//     (
//         useTimeName_
//       ? time_.timeName()
//       : printf_.empty()
//       ? Foam::name(time_.timeIndex())
//       : word::printf(printf_, time_.timeIndex())
//     );

    const word timeDesc =
    (
        printf_.empty()
      ? Foam::name(time_.timeIndex())
      : word::printf(printf_, time_.timeIndex())
    );

    fileName vtkDir(dirName_);
    vtkDir.expand();
    if (!vtkDir.isAbsolute())
    {
        vtkDir = stringOps::expand("<case>")/vtkDir;
    }
    mkDir(vtkDir);

    Log << name() << " output Time: " << time_.timeName() << nl;

    // Each cloud separately
    for (const word& cloudName : cloudNames)
    {
        // Legacy is not to be supported

        const fileName outputName
        (
            vtkDir/cloudName + "_" + timeDesc + ".vtp"
        );

        if (writeCloud(outputName, cloudName))
        {
            Log << "    cloud  : " << outputName << endl;

            if (Pstream::master())
            {
                // Add to file-series and emit as JSON
                // - misbehaves if vtkDir changes during the run,
                // but that causes other issues too.

                series_(cloudName).append({time_.value(), timeDesc});

                vtk::seriesWrite
                (
                    vtkDir/cloudName + ".vtp",
                    series_[cloudName]
                );
            }
        }
    }

    return true;
}


// ************************************************************************* //
