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

Application
    foamToVTK

Group
    grpPostProcessingUtilities

Description
    VTK file format writer.

    - Handles volFields, pointFields, surfaceScalarField, surfaceVectorField
      fields.
    - Mesh topo changes.
    - Both ascii and binary.
    - Single time step writing.
    - Write subset only.
    - Automatic decomposition of cells; polygons on boundary undecomposed since
      handled by vtk.

Usage
    \b foamToVTK [OPTION]

    Options:
      - \par -ascii
        Write VTK data in ASCII format instead of binary.

      - \par -xml
        Write VTK data in XML format instead of legacy format

      - \par -mesh \<name\>
        Use a different mesh name (instead of -region)

      - \par -fields \<fields\>
        Convert selected fields only. For example,
        \verbatim
          -fields '( p T U )'
        \endverbatim
        The quoting is required to avoid shell expansions and to pass the
        information as a single argument.

      - \par -surfaceFields
        Write surfaceScalarFields (e.g., phi)

      - \par -cellSet \<name\>
      - \par -cellZone \<name\>
        Restrict conversion to either the cellSet or the cellZone.

      - \par -faceSet \<name\>
      - \par -pointSet \<name\>
        Restrict conversion to the faceSet or pointSet.

      - \par -nearCellValue
        Output cell value on patches instead of patch value itself

      - \par -noInternal
        Do not generate file for mesh, only for patches

      - \par -noLagrangian
        Suppress writing lagrangian positions and fields.

      - \par -noPointValues
        No pointFields

      - \par -noFaceZones
        No faceZones

      - \par -noLinks
        (in parallel) do not link processor files to master

      - \par poly
        write polyhedral cells without tet/pyramid decomposition

      - \par -allPatches
        Combine all patches into a single file

      - \par -excludePatches \<patchNames\>
        Specify patches (wildcards) to exclude. For example,
        \verbatim
          -excludePatches '( inlet_1 inlet_2 "proc.*")'
        \endverbatim
        The quoting is required to avoid shell expansions and to pass the
        information as a single argument. The double quotes denote a regular
        expression.

      - \par -useTimeName
        use the time index in the VTK file name instead of the time index

Note
    mesh subset is handled by meshSubsetHelper. Slight inconsistency in
    interpolation: on the internal field it interpolates the whole volField
    to the whole-mesh pointField and then selects only those values it
    needs for the subMesh (using the fvMeshSubset cellMap(), pointMap()
    functions). For the patches however it uses the
    fvMeshSubset.interpolate function to directly interpolate the
    whole-mesh values onto the subset patch.

Note
    \par new file format:
    no automatic timestep recognition.
    However can have .pvd file format which refers to time simulation
    if XML *.vtu files are available:

    \verbatim
      <?xml version="1.0"?>
      <VTKFile type="Collection" version="0.1" byte_order="LittleEndian">
        <Collection>
          <DataSet timestep="50" file="pitzDaily_2.vtu"/>
          <DataSet timestep="100" file="pitzDaily_3.vtu"/>
          <DataSet timestep="150" file="pitzDaily_4.vtu"/>
          <DataSet timestep="200" file="pitzDaily_5.vtu"/>
          <DataSet timestep="250" file="pitzDaily_6.vtu"/>
          <DataSet timestep="300" file="pitzDaily_7.vtu"/>
          <DataSet timestep="350" file="pitzDaily_8.vtu"/>
          <DataSet timestep="400" file="pitzDaily_9.vtu"/>
          <DataSet timestep="450" file="pitzDaily_10.vtu"/>
          <DataSet timestep="500" file="pitzDaily_11.vtu"/>
        </Collection>
      </VTKFile>
    \endverbatim

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pointMesh.H"
#include "volPointInterpolation.H"
#include "emptyPolyPatch.H"
#include "labelIOField.H"
#include "scalarIOField.H"
#include "sphericalTensorIOField.H"
#include "symmTensorIOField.H"
#include "tensorIOField.H"
#include "faceZoneMesh.H"
#include "Cloud.H"
#include "passiveParticle.H"
#include "stringOps.H"
#include "areaFields.H"

#include "meshSubsetHelper.H"
#include "readFields.H"
#include "faceSet.H"
#include "pointSet.H"

#include "foamVtkOutputOptions.H"
#include "foamVtkInternalWriter.H"
#include "foamVtkPatchWriter.H"
#include "foamVtkSurfaceMeshWriter.H"
#include "foamVtkLagrangianWriter.H"
#include "foamVtkWriteFaceSet.H"
#include "foamVtkWritePointSet.H"
#include "foamVtkWriteSurfFields.H"

#include "memInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoField>
void print(const char* msg, Ostream& os, const UPtrList<const GeoField>& flds)
{
    if (flds.size())
    {
        os  << msg;
        forAll(flds, i)
        {
            os  << ' ' << flds[i].name();
        }
        os  << endl;
    }
}


void print(Ostream& os, const wordList& flds)
{
    forAll(flds, i)
    {
        os  << ' ' << flds[i];
    }
    os  << endl;
}


labelList getSelectedPatches
(
    const polyBoundaryMesh& patches,
    const wordRes& excludePatches
)
{
    DynamicList<label> patchIDs(patches.size());

    Info<< "Combining patches:" << endl;

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if
        (
            isType<emptyPolyPatch>(pp)
         || (Pstream::parRun() && isType<processorPolyPatch>(pp))
        )
        {
            Info<< "    discarding empty/processor patch " << patchi
                << " " << pp.name() << endl;
        }
        else if (excludePatches.match(pp.name()))
        {
            Info<< "    excluding patch " << patchi
                << " " << pp.name() << endl;
        }
        else
        {
            patchIDs.append(patchi);
            Info<< "    patch " << patchi << " " << pp.name() << endl;
        }
    }

    return patchIDs.shrink();
}


HashTable<wordHashSet> candidateObjects
(
    const IOobjectList& objects,
    const wordHashSet& supportedTypes,
    const bool specifiedFields,
    const wordHashSet& selectedFields
)
{
    // Special case = no fields
    if (specifiedFields && selectedFields.empty())
    {
        return HashTable<wordHashSet>();
    }

    HashTable<wordHashSet> usable = objects.classes();

    // Limited to types that we explicitly handle
    usable.retain(supportedTypes);

    // If specified, further limit to selected fields
    if (specifiedFields)
    {
        forAllIters(usable, iter)
        {
            iter.object().retain(selectedFields);
        }

        // Prune entries without any fields
        usable.filterValues
        (
            [](const wordHashSet& vals){ return !vals.empty(); }
        );
    }

    return usable;
}


//
// Process args for output options
// Default from foamVtkOutputOptions is inline ASCII xml
//
vtk::outputOptions getOutputOptions(const argList& args)
{
    vtk::outputOptions opts;

    if (args.optionFound("xml"))
    {
        opts.ascii(args.optionFound("ascii"));
    }
    else
    {
        opts.legacy(true);

        if (!args.optionFound("ascii"))
        {
            if (sizeof(floatScalar) != 4 || sizeof(label) != 4)
            {
                opts.ascii(true);

                WarningInFunction
                    << "Using ASCII rather than legacy binary VTK format since "
                    << "floatScalar and/or label are not 4 bytes in size."
                    << nl << endl;
            }
            else
            {
                opts.ascii(false);
            }
        }
    }

    return opts;
}


fileName relativeName(const fileName& parent, const fileName& file)
{
    string::size_type next = parent.size();
    if
    (
        file.startsWith(parent)
     && next < file.size()
     && file[next] == '/'
    )
    {
        return file.substr(next+1);
    }
    else
    {
        return file;
    }
}


fileName relativeName(const Time& runTime, const fileName& file)
{
    return relativeName(runTime.path(), file);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote("VTK file format writer");
    timeSelector::addOptions();

    #include "addRegionOption.H"

    argList::addOption
    (
        "fields",
        "wordList",
        "only convert the specified fields - eg '(p T U)'"
    );
    argList::addOption
    (
        "cellSet",
        "name",
        "convert a mesh subset corresponding to the specified cellSet"
    );
    argList::addOption
    (
        "cellZone",
        "name",
        "convert a mesh subset corresponding to the specified cellZone"
    );
    argList::addOption
    (
        "faceSet",
        "name",
        "restrict conversion to the specified faceSet"
    );
    argList::addOption
    (
        "pointSet",
        "name",
        "restrict conversion to the specified pointSet"
    );
    argList::addBoolOption
    (
        "ascii",
        "write in ASCII format instead of binary"
    );
    argList::addBoolOption
    (
        "xml",
        "write VTK xml instead of legacy format"
    );
    argList::addBoolOption
    (
        "poly",
        "write polyhedral cells without tet/pyramid decomposition"
    );
    argList::addBoolOption
    (
        "surfaceFields",
        "write surfaceScalarFields (e.g., phi)"
    );
    argList::addBoolOption
    (
        "finiteAreaFields",
        "write finite area fields"
    );
    argList::addBoolOption
    (
        "nearCellValue",
        "use cell value on patches instead of patch value itself"
    );
    argList::addBoolOption
    (
        "noInternal",
        "do not generate file for mesh, only for patches"
    );
    argList::addBoolOption
    (
        "noLagrangian",
        "suppress writing lagrangian positions and fields"
    );
    argList::addBoolOption
    (
        "noPointValues",
        "no pointFields"
    );
    argList::addBoolOption
    (
        "allPatches",
        "combine all patches into a single file"
    );
    argList::addOption
    (
        "excludePatches",
        "wordReList",
        "a list of patches to exclude - eg '( inlet \".*Wall\" )' "
    );
    argList::addBoolOption
    (
        "noFaceZones",
        "no faceZones"
    );
    argList::addBoolOption
    (
        "noLinks",
        "don't link processor VTK files - parallel only"
    );
    argList::addBoolOption
    (
        "useTimeName",
        "use time name instead of the time index when naming files"
    );
    argList::addOption
    (
        "name",
        "subdir",
        "sub-directory name for VTK output (default: 'VTK')"
    );

    #include "setRootCase.H"

    cpuTime timer;
    memInfo mem;
    Info<< "Initial memory "
        << mem.update().size() << " kB" << endl;

    #include "createTime.H"

    const bool decomposePoly   = !args.optionFound("poly");
    const bool doWriteInternal = !args.optionFound("noInternal");
    const bool doFaceZones     = !args.optionFound("noFaceZones");
    const bool doLinks         = !args.optionFound("noLinks");
    const bool useTimeName     = args.optionFound("useTimeName");
    const bool noLagrangian    = args.optionFound("noLagrangian");
    const bool nearCellValue   = args.optionFound("nearCellValue");

    const vtk::outputOptions fmtType = getOutputOptions(args);

    if (nearCellValue)
    {
        WarningInFunction
            << "Using neighbouring cell value instead of patch value"
            << nl << endl;
    }

    const bool noPointValues = args.optionFound("noPointValues");

    if (noPointValues)
    {
        Info<< "Outputting cell values only."
            << " Point fields disabled by '-noPointValues' option"
            << nl;
    }

    const bool allPatches = args.optionFound("allPatches");

    wordReList excludePatches;
    if (args.optionFound("excludePatches"))
    {
        args.optionLookup("excludePatches")() >> excludePatches;

        Info<< "Not including patches " << excludePatches << nl << endl;
    }

    string vtkName = runTime.caseName();

    meshSubsetHelper::subsetType cellSubsetType = meshSubsetHelper::NONE;
    word cellSubsetName;
    if (args.optionReadIfPresent("cellSet", cellSubsetName))
    {
        vtkName = cellSubsetName;
        cellSubsetType = meshSubsetHelper::SET;
    }
    else if (args.optionReadIfPresent("cellZone", cellSubsetName))
    {
        vtkName = cellSubsetName;
        cellSubsetType = meshSubsetHelper::ZONE;
    }
    else if (Pstream::parRun())
    {
        // Strip off leading casename, leaving just processor_DDD ending.
        string::size_type i = vtkName.rfind("processor");

        if (i != string::npos)
        {
            vtkName = vtkName.substr(i);
        }
    }

    word faceSetName;
    args.optionReadIfPresent("faceSet", faceSetName);

    word pointSetName;
    args.optionReadIfPresent("pointSet", pointSetName);

    // Define sub-directory name to use for VTK data.
    const word vtkDirName = args.optionLookupOrDefault<word>("name", "VTK");

    #include "createNamedMesh.H"

    // VTK/ directory in the case
    fileName fvPath(runTime.path()/vtkDirName);

    // Directory of mesh (region0 gets filtered out)
    fileName regionPrefix;
    if (regionName != polyMesh::defaultRegion)
    {
        fvPath = fvPath/regionName;
        regionPrefix = regionName;
    }

    if (isDir(fvPath))
    {
        if
        (
            args.optionFound("time")
         || args.optionFound("latestTime")
         || cellSubsetName.size()
         || faceSetName.size()
         || pointSetName.size()
         || regionName != polyMesh::defaultRegion
        )
        {
            Info<< "Keeping old VTK files in " << fvPath << nl << endl;
        }
        else
        {
            Info<< "Deleting old VTK files in " << fvPath << nl << endl;

            rmDir(fvPath);
        }
    }

    mkDir(fvPath);

    instantList timeDirs = timeSelector::select0(runTime, args);

    // Mesh wrapper: does subsetting and decomposition
    meshSubsetHelper meshRef(mesh, cellSubsetType, cellSubsetName);

    // Collect decomposition information etc.
    vtk::vtuCells vtuMeshCells(fmtType, decomposePoly);

    Info<< "VTK mesh topology: "
        << timer.cpuTimeIncrement() << " s, "
        << mem.update().size() << " kB" << endl;

    #include "findClouds.H"

    // Supported volume field types
    const wordHashSet vFieldTypes
    {
        volScalarField::typeName,
        volVectorField::typeName,
        volSphericalTensorField::typeName,
        volSymmTensorField::typeName,
        volTensorField::typeName
    };

    // Supported dimensioned field types
    const wordHashSet dFieldTypes
    {
        volScalarField::Internal::typeName,
        volVectorField::Internal::typeName,
        volSphericalTensorField::Internal::typeName,
        volSymmTensorField::Internal::typeName,
        volTensorField::Internal::typeName
    };

    // Supported point field types
    const wordHashSet pFieldTypes
    {
        pointScalarField::typeName,
        pointVectorField::typeName,
        pointSphericalTensorField::typeName,
        pointSymmTensorField::typeName,
        pointTensorField::typeName
    };

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);

        Info<< "Time: " << runTime.timeName() << endl;

        const word timeDesc =
            useTimeName ? runTime.timeName() : Foam::name(runTime.timeIndex());

        // Check for new polyMesh/ and update mesh, fvMeshSubset and cell
        // decomposition.
        polyMesh::readUpdateState meshState = meshRef.readUpdate();

        const fvMesh& mesh = meshRef.mesh();
        if
        (
            meshState == polyMesh::TOPO_CHANGE
         || meshState == polyMesh::TOPO_PATCH_CHANGE
        )
        {
            // Trigger change for vtk cells too
            vtuMeshCells.clear();
        }

        // If faceSet: write faceSet only (as polydata)
        if (faceSetName.size())
        {
            // Load the faceSet
            faceSet set(mesh, faceSetName);

            // Filename as if patch with same name.
            mkDir(fvPath/set.name());

            fileName outputName
            (
                fvPath/set.name()/set.name()
              + "_"
              + timeDesc
            );
            Info<< "    faceSet   : "
                << relativeName(runTime, outputName) << nl;

            vtk::writeFaceSet
            (
                meshRef.mesh(),
                set,
                outputName,
                fmtType
            );
            continue;
        }

        // If pointSet: write pointSet only (as polydata)
        if (pointSetName.size())
        {
            // Load the pointSet
            pointSet set(mesh, pointSetName);

            // Filename as if patch with same name.
            mkDir(fvPath/set.name());

            fileName outputName
            (
                fvPath/set.name()/set.name()
              + "_"
              + timeDesc
            );
            Info<< "    pointSet  : "
                << relativeName(runTime, outputName) << nl;

            vtk::writePointSet
            (
                meshRef.mesh(),
                set,
                outputName,
                fmtType
            );
            continue;
        }


        // Search for list of objects for this time
        IOobjectList objects(mesh, runTime.timeName());

        wordHashSet selectedFields;
        const bool specifiedFields = args.optionReadIfPresent
        (
            "fields",
            selectedFields
        );

        // Construct the vol fields
        // References the original mesh, but uses subsetted portion only.

        PtrList<const volScalarField> vScalarFld;
        PtrList<const volVectorField> vVectorFld;
        PtrList<const volSphericalTensorField> vSphTensorf;
        PtrList<const volSymmTensorField> vSymTensorFld;
        PtrList<const volTensorField> vTensorFld;

        if
        (
            candidateObjects
            (
                objects,
                vFieldTypes,
                specifiedFields,
                selectedFields
            ).size()
        )
        {
            readFields
            (
                meshRef,
                meshRef.baseMesh(),
                objects,
                selectedFields,
                vScalarFld
            );
            print("    volScalar        :", Info, vScalarFld);

            readFields
            (
                meshRef,
                meshRef.baseMesh(),
                objects,
                selectedFields,
                vVectorFld
            );
            print("    volVector        :", Info, vVectorFld);

            readFields
            (
                meshRef,
                meshRef.baseMesh(),
                objects,
                selectedFields,
                vSphTensorf
            );
            print("    volSphTensor     :", Info, vSphTensorf);

            readFields
            (
                meshRef,
                meshRef.baseMesh(),
                objects,
                selectedFields,
                vSymTensorFld
            );
            print("    volSymmTensor    :", Info, vSymTensorFld);

            readFields
            (
                meshRef,
                meshRef.baseMesh(),
                objects,
                selectedFields,
                vTensorFld
            );
            print("    volTensor        :", Info, vTensorFld);
        }

        const label nVolFields =
        (
            vScalarFld.size()
          + vVectorFld.size()
          + vSphTensorf.size()
          + vSymTensorFld.size()
          + vTensorFld.size()
        );


        // Construct dimensioned fields
        PtrList<const volScalarField::Internal> dScalarFld;
        PtrList<const volVectorField::Internal> dVectorFld;
        PtrList<const volSphericalTensorField::Internal> dSphTensorFld;
        PtrList<const volSymmTensorField::Internal> dSymTensorFld;
        PtrList<const volTensorField::Internal> dTensorFld;

        if
        (
            candidateObjects
            (
                objects,
                dFieldTypes,
                specifiedFields,
                selectedFields
            ).size()
        )
        {
            readFields
            (
                meshRef,
                meshRef.baseMesh(),
                objects,
                selectedFields,
                dScalarFld
            );
            print("    volScalar::Internal      :", Info, dScalarFld);

            readFields
            (
                meshRef,
                meshRef.baseMesh(),
                objects,
                selectedFields,
                dVectorFld
            );
            print("    volVector::Internal      :", Info, dVectorFld);

            readFields
            (
                meshRef,
                meshRef.baseMesh(),
                objects,
                selectedFields,
                dSphTensorFld
            );
            print("    volSphTensor::Internal   :", Info, dSphTensorFld);

            readFields
            (
                meshRef,
                meshRef.baseMesh(),
                objects,
                selectedFields,
                dSymTensorFld
            );
            print("    volSymmTensor::Internal  :", Info, dSymTensorFld);

            readFields
            (
                meshRef,
                meshRef.baseMesh(),
                objects,
                selectedFields,
                dTensorFld
            );
            print("    volTensor::Internal      :", Info, dTensorFld);
        }

        const label nDimFields =
        (
            dScalarFld.size()
          + dVectorFld.size()
          + dSphTensorFld.size()
          + dSymTensorFld.size()
          + dTensorFld.size()
        );


        // Finite-area mesh and fields - need not exist

        if (args.optionFound("finiteAreaFields"))
        {
            autoPtr<faMesh> aMeshPtr;
            {
                const bool throwing = FatalError.throwExceptions();
                try
                {
                    aMeshPtr.reset(new faMesh(meshRef.baseMesh()));
                }
                catch (Foam::error& err)
                {
                    aMeshPtr.clear();
                }
                FatalError.throwExceptions(throwing);
            }

            if (aMeshPtr.valid())
            {
                // Construct the area fields

                PtrList<const areaScalarField> aScalarFld;
                PtrList<const areaVectorField> aVectorFld;
                PtrList<const areaSphericalTensorField> aSphTensorf;
                PtrList<const areaSymmTensorField> aSymTensorFld;
                PtrList<const areaTensorField> aTensorFld;

                const faMesh& aMesh = aMeshPtr();

                if (!specifiedFields || selectedFields.size())
                {
                    readFields
                    (
                        aMesh,
                        objects,
                        selectedFields,
                        aScalarFld
                    );
                    print("    areaScalar           :", Info, aScalarFld);

                    readFields
                    (
                        aMesh,
                        objects,
                        selectedFields,
                        aVectorFld
                    );
                    print("    areaVector           :", Info, aVectorFld);

                    readFields
                    (
                        aMesh,
                        objects,
                        selectedFields,
                        aSphTensorf
                    );
                    print("    areaSphericalTensor   :", Info, aSphTensorf);

                    readFields
                    (
                        aMesh,
                        objects,
                        selectedFields,
                        aSymTensorFld
                    );
                    print("    areaSymmTensor        :", Info, aSymTensorFld);

                    readFields
                    (
                        aMesh,
                        objects,
                        selectedFields,
                        aTensorFld
                    );
                    print("    areaTensor            :", Info, aTensorFld);
                }

                const label nAreaFields =
                (
                    aScalarFld.size()
                  + aVectorFld.size()
                  + aSphTensorf.size()
                  + aSymTensorFld.size()
                  + aTensorFld.size()
                );

                fileName outputName(fvPath/"finiteArea");

                mkDir(outputName);

                const auto& pp = aMesh.patch();

                vtk::surfaceMeshWriter writer
                (
                    pp,
                    aMesh.name(),
                    outputName/"finiteArea" + "_" + timeDesc,
                    fmtType
                );

                // Number of fields
                writer.beginCellData(nAreaFields);

                writer.write(aScalarFld);
                writer.write(aVectorFld);
                writer.write(aSphTensorf);
                writer.write(aSymTensorFld);
                writer.write(aTensorFld);

                writer.endCellData();

                writer.writeFooter();
            }
        }

        PtrList<const pointScalarField> pScalarFld;
        PtrList<const pointVectorField> pVectorFld;
        PtrList<const pointSphericalTensorField> pSphTensorFld;
        PtrList<const pointSymmTensorField> pSymTensorFld;
        PtrList<const pointTensorField> pTensorFld;

        // Construct pointMesh only if necessary since it constructs edge
        // addressing (expensive on polyhedral meshes)
        if
        (
            !noPointValues
         && candidateObjects
            (
                objects,
                pFieldTypes,
                specifiedFields,
                selectedFields
            ).size()
        )
        {
            const pointMesh& ptMesh = pointMesh::New(meshRef.baseMesh());

            readFields
            (
                meshRef,
                ptMesh,
                objects,
                selectedFields,
                pScalarFld
            );
            print("    pointScalar      :", Info, pScalarFld);

            readFields
            (
                meshRef,
                ptMesh,
                objects,
                selectedFields,
                pVectorFld
            );
            print("    pointVector      :", Info, pVectorFld);

            readFields
            (
                meshRef,
                ptMesh,
                objects,
                selectedFields,
                pSphTensorFld
            );
            print("    pointSphTensor   : ", Info, pSphTensorFld);

            readFields
            (
                meshRef,
                ptMesh,
                objects,
                selectedFields,
                pSymTensorFld
            );
            print("    pointSymmTensor  :", Info, pSymTensorFld);

            readFields
            (
                meshRef,
                ptMesh,
                objects,
                selectedFields,
                pTensorFld
            );
            print("    pointTensor      :", Info, pTensorFld);
        }

        const label nPointFields =
            pScalarFld.size()
          + pVectorFld.size()
          + pSphTensorFld.size()
          + pSymTensorFld.size()
          + pTensorFld.size();

        if (doWriteInternal)
        {
            if (vtuMeshCells.empty())
            {
                // subMesh or baseMesh
                vtuMeshCells.reset(meshRef.mesh());
            }

            // Create file and write header
            fileName outputName
            (
                fvPath/vtkName
              + "_"
              + timeDesc
            );
            Info<< "    Internal  : "
                << relativeName(runTime, outputName) << endl;

            // Write mesh
            vtk::internalWriter writer
            (
                meshRef.mesh(),
                vtuMeshCells,
                outputName,
                fmtType
            );

            // CellData
            {
                writer.beginCellData(1 + nVolFields + nDimFields);

                // Write cellID field
                writer.writeCellIDs();

                // Write volFields
                writer.write(vScalarFld);
                writer.write(vVectorFld);
                writer.write(vSphTensorf);
                writer.write(vSymTensorFld);
                writer.write(vTensorFld);

                // Write dimensionedFields
                writer.write(dScalarFld);
                writer.write(dVectorFld);
                writer.write(dSphTensorFld);
                writer.write(dSymTensorFld);
                writer.write(dTensorFld);

                writer.endCellData();
            }

            // PointData
            if (!noPointValues)
            {
                writer.beginPointData(nVolFields + nDimFields + nPointFields);

                // pointFields
                writer.write(pScalarFld);
                writer.write(pVectorFld);
                writer.write(pSphTensorFld);
                writer.write(pSymTensorFld);
                writer.write(pTensorFld);

                // Interpolated volFields
                volPointInterpolation pInterp(mesh);

                writer.write(pInterp, vScalarFld);
                writer.write(pInterp, vVectorFld);
                writer.write(pInterp, vSphTensorf);
                writer.write(pInterp, vSymTensorFld);
                writer.write(pInterp, vTensorFld);

                writer.write(pInterp, dScalarFld);
                writer.write(pInterp, dVectorFld);
                writer.write(pInterp, dSphTensorFld);
                writer.write(pInterp, dSymTensorFld);
                writer.write(pInterp, dTensorFld);

                writer.endPointData();
            }

            writer.writeFooter();
        }

        //---------------------------------------------------------------------
        //
        // Write surface fields
        //
        //---------------------------------------------------------------------

        if (args.optionFound("surfaceFields"))
        {
            PtrList<const surfaceScalarField> sScalarFld;
            readFields
            (
                meshRef,
                meshRef.baseMesh(),
                objects,
                selectedFields,
                sScalarFld
            );
            print("    surfScalar   :", Info, sScalarFld);

            PtrList<const surfaceVectorField> sVectorFld;
            readFields
            (
                meshRef,
                meshRef.baseMesh(),
                objects,
                selectedFields,
                sVectorFld
            );
            print("    surfVector   :", Info, sVectorFld);

            if (sScalarFld.size())
            {
                // Rework the scalar fields into vector fields.
                const label sz = sVectorFld.size();

                sVectorFld.setSize(sz + sScalarFld.size());

                surfaceVectorField n(mesh.Sf()/mesh.magSf());

                forAll(sScalarFld, i)
                {
                    surfaceVectorField* ssfPtr = (sScalarFld[i]*n).ptr();
                    ssfPtr->rename(sScalarFld[i].name());
                    sVectorFld.set(sz+i, ssfPtr);
                }
                sScalarFld.clear();
            }

            if (sVectorFld.size())
            {
                mkDir(fvPath / "surfaceFields");

                fileName outputName
                (
                    fvPath
                  / "surfaceFields"
                  / "surfaceFields"
                  + "_"
                  + timeDesc
                );

                vtk::writeSurfFields
                (
                    meshRef.mesh(),
                    outputName,
                    fmtType,
                    sVectorFld
                );
            }
        }


        //---------------------------------------------------------------------
        //
        // Write patches (POLYDATA file, one for each patch)
        //
        //---------------------------------------------------------------------

        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        if (allPatches)
        {
            mkDir(fvPath/"allPatches");

            fileName outputName
            (
                fvPath/"allPatches"
              / (meshRef.useSubMesh() ? cellSubsetName : "allPatches")
              + "_"
              + timeDesc
            );
            Info<< "    Combined patches    : "
                << relativeName(runTime, outputName) << nl;

            vtk::patchWriter writer
            (
                meshRef.mesh(),
                outputName,
                fmtType,
                nearCellValue,
                getSelectedPatches(patches, excludePatches)
            );

            // CellData
            {
                writer.beginCellData(1 + nVolFields);

                // Write patchID field
                writer.writePatchIDs();

                // Write volFields
                writer.write(vScalarFld);
                writer.write(vVectorFld);
                writer.write(vSphTensorf);
                writer.write(vSymTensorFld);
                writer.write(vTensorFld);

                writer.endCellData();
            }

            // PointData
            if (!noPointValues)
            {
                writer.beginPointData(nPointFields);

                // Write pointFields
                writer.write(pScalarFld);
                writer.write(pVectorFld);
                writer.write(pSphTensorFld);
                writer.write(pSymTensorFld);
                writer.write(pTensorFld);

                // no interpolated volFields to avoid creating
                // patchInterpolation for all subpatches.

                writer.endPointData();
            }

            writer.writeFooter();
        }
        else
        {
            forAll(patches, patchi)
            {
                const polyPatch& pp = patches[patchi];

                if (stringOps::match(excludePatches, pp.name()))
                {
                    // Skip excluded patch
                    continue;
                }

                mkDir(fvPath/pp.name());

                fileName outputName
                (
                    fvPath/pp.name()
                  / (meshRef.useSubMesh() ? cellSubsetName : pp.name())
                  + "_"
                  + timeDesc
                );
                Info<< "    Patch     : "
                    << relativeName(runTime, outputName) << nl;

                vtk::patchWriter writer
                (
                    meshRef.mesh(),
                    outputName,
                    fmtType,
                    nearCellValue,
                    labelList{patchi}
                );

                if (!isA<emptyPolyPatch>(pp))
                {
                    // VolFields + patchID
                    writer.beginCellData(1 + nVolFields);

                    // Write patchID field
                    writer.writePatchIDs();

                    // Write volFields
                    writer.write(vScalarFld);
                    writer.write(vVectorFld);
                    writer.write(vSphTensorf);
                    writer.write(vSymTensorFld);
                    writer.write(vTensorFld);

                    writer.endCellData();

                    if (!noPointValues)
                    {
                        writer.beginPointData(nVolFields + nPointFields);

                        // Write pointFields
                        writer.write(pScalarFld);
                        writer.write(pVectorFld);
                        writer.write(pSphTensorFld);
                        writer.write(pSymTensorFld);
                        writer.write(pTensorFld);

                        PrimitivePatchInterpolation<primitivePatch> pInter
                        (
                            pp
                        );

                        // Write interpolated volFields
                        writer.write(pInter, vScalarFld);
                        writer.write(pInter, vVectorFld);
                        writer.write(pInter, vSphTensorf);
                        writer.write(pInter, vSymTensorFld);
                        writer.write(pInter, vTensorFld);

                        writer.endPointData();
                    }
                }

                writer.writeFooter();
            }
        }

        //---------------------------------------------------------------------
        //
        // Write faceZones (POLYDATA file, one for each zone)
        //
        //---------------------------------------------------------------------

        if (doFaceZones && !mesh.faceZones().empty())
        {
            PtrList<const surfaceScalarField> sScalarFld;
            readFields
            (
                meshRef,
                meshRef.baseMesh(),
                objects,
                selectedFields,
                sScalarFld
            );
            print("    surfScalar   :", Info, sScalarFld);

            PtrList<const surfaceVectorField> sVectorFld;
            readFields
            (
                meshRef,
                meshRef.baseMesh(),
                objects,
                selectedFields,
                sVectorFld
            );
            print("    surfVector   :", Info, sVectorFld);

            for (const faceZone& fz : mesh.faceZones())
            {
                mkDir(fvPath/fz.name());

                fileName outputName =
                (
                    fvPath/fz.name()
                  / (meshRef.useSubMesh() ? cellSubsetName : fz.name())
                  + "_"
                  + timeDesc
                );
                Info<< "    FaceZone  : "
                    << relativeName(runTime, outputName) << nl;

                indirectPrimitivePatch pp
                (
                    IndirectList<face>(mesh.faces(), fz),
                    mesh.points()
                );

                vtk::surfaceMeshWriter writer
                (
                    pp,
                    fz.name(),
                    outputName,
                    fmtType
                );

                // Number of fields
                writer.beginCellData(sScalarFld.size() + sVectorFld.size());

                writer.write(sScalarFld);
                writer.write(sVectorFld);

                writer.endCellData();

                writer.writeFooter();
            }
        }


        //---------------------------------------------------------------------
        //
        // Write lagrangian data
        //
        //---------------------------------------------------------------------

        for (const fileName& cloudName : cloudNames)
        {
            // Always create the cloud directory.
            mkDir(fvPath/cloud::prefix/cloudName);

            fileName outputName
            (
                fvPath/cloud::prefix/cloudName/cloudName
              + "_" + timeDesc
            );
            Info<< "    Lagrangian: "
                << relativeName(runTime, outputName) << nl;

            IOobjectList sprayObjs
            (
                mesh,
                runTime.timeName(),
                cloud::prefix/cloudName
            );

            if (sprayObjs.found("positions") || sprayObjs.found("coordinates"))
            {
                wordList labelNames(sprayObjs.names(labelIOField::typeName));
                Info<< "        labels      :";
                print(Info, labelNames);

                wordList scalarNames(sprayObjs.names(scalarIOField::typeName));
                Info<< "        scalars     :";
                print(Info, scalarNames);

                wordList vectorNames(sprayObjs.names(vectorIOField::typeName));
                Info<< "        vectors     :";
                print(Info, vectorNames);

                wordList sphereNames
                (
                    sprayObjs.names
                    (
                        sphericalTensorIOField::typeName
                    )
                );
                Info<< "        sphTensors  :";
                print(Info, sphereNames);

                wordList symmNames
                (
                    sprayObjs.names
                    (
                        symmTensorIOField::typeName
                    )
                );
                Info<< "        symmTensors :";
                print(Info, symmNames);

                wordList tensorNames(sprayObjs.names(tensorIOField::typeName));
                Info<< "        tensors     :";
                print(Info, tensorNames);

                vtk::lagrangianWriter writer
                (
                    meshRef.mesh(),
                    cloudName,
                    outputName,
                    fmtType
                );

                // Write number of fields
                writer.beginParcelData
                (
                    labelNames.size()
                  + scalarNames.size()
                  + vectorNames.size()
                  + sphereNames.size()
                  + symmNames.size()
                  + tensorNames.size()
                );

                // Fields
                writer.writeIOField<label>(labelNames);
                writer.writeIOField<scalar>(scalarNames);
                writer.writeIOField<vector>(vectorNames);
                writer.writeIOField<sphericalTensor>(sphereNames);
                writer.writeIOField<symmTensor>(symmNames);
                writer.writeIOField<tensor>(tensorNames);

                writer.endParcelData();

                writer.writeFooter();
            }
            else
            {
                vtk::lagrangianWriter writer
                (
                    meshRef.mesh(),
                    cloudName,
                    outputName,
                    fmtType,
                    true
                );

                // Write number of fields
                writer.beginParcelData(0);

                writer.endParcelData();

                writer.writeFooter();
            }
        }

        Info<< "Wrote in "
            << timer.cpuTimeIncrement() << " s, "
            << mem.update().size() << " kB" << endl;
    }


    //---------------------------------------------------------------------
    //
    // Link parallel outputs back to undecomposed case for ease of loading
    //
    //---------------------------------------------------------------------

    if (Pstream::parRun() && doLinks)
    {
        mkDir(runTime.path()/".."/vtkDirName);
        chDir(runTime.path()/".."/vtkDirName);

        Info<< "Linking all processor files to "
            << runTime.path()/".."/vtkDirName
            << endl;

        // Get list of vtk files
        fileName procVTK
        (
            fileName("..")
          / "processor" + Foam::name(Pstream::myProcNo())
          / vtkDirName
        );

        fileNameList dirs(readDir(procVTK, fileName::DIRECTORY));
        label sz = dirs.size();
        dirs.setSize(sz+1);
        dirs[sz] = ".";

        for (const fileName& subDir : dirs)
        {
            fileNameList subFiles(readDir(procVTK/subDir, fileName::FILE));

            for (const fileName& subFile : subFiles)
            {
                fileName procFile(procVTK/subDir/subFile);

                if (exists(procFile))
                {
                    // Could likely also use Foam::ln() directly
                    List<string> cmd
                    {
                        "ln",
                        "-s",
                        procFile,
                        (
                            "processor"
                          + Foam::name(Pstream::myProcNo())
                          + "_"
                          + procFile.name()
                        )
                    };

                    if (Foam::system(cmd) == -1)
                    {
                        WarningInFunction
                            << "Could not execute command " << cmd << endl;
                    }
                }
            }
        }
    }

    Info<< "\nEnd: "
        << timer.elapsedCpuTime() << " s, "
        << mem.update().peak() << " kB (peak)\n" << endl;

    return 0;
}


// ************************************************************************* //
