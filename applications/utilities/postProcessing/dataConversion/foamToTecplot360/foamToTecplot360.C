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

Application
    foamToTecplot360

Group
    grpPostProcessingUtilities

Description
    Tecplot binary file format writer.

Usage
    \b foamToTecplot360 [OPTION]

    Options:
      - \par -fields \<names\>
        Convert selected fields only. For example,
        \verbatim
          -fields '( p T U )'
        \endverbatim
        The quoting is required to avoid shell expansions and to pass the
        information as a single argument.

      - \par -cellSet \<name\>
      - \par -faceSet \<name\>
        Restrict conversion to the cellSet, faceSet.

      - \par -nearCellValue
        Output cell value on patches instead of patch value itself

      - \par -noInternal
        Do not generate file for mesh, only for patches

      - \par -noPointValues
        No pointFields

      - \par -noFaceZones
        No faceZones

      - \par -excludePatches \<patchNames\>
        Specify patches (wildcards) to exclude. For example,
        \verbatim
          -excludePatches '( inlet_1 inlet_2 "proc.*")'
        \endverbatim
        The quoting is required to avoid shell expansions and to pass the
        information as a single argument. The double quotes denote a regular
        expression.

\*---------------------------------------------------------------------------*/

#include "pointMesh.H"
#include "volPointInterpolation.H"
#include "emptyPolyPatch.H"
#include "labelIOField.H"
#include "scalarIOField.H"
#include "sphericalTensorIOField.H"
#include "symmTensorIOField.H"
#include "tensorIOField.H"
#include "passiveParticleCloud.H"
#include "faceSet.H"
#include "stringOps.H"
#include "wordReList.H"

#include "meshSubsetHelper.H"
#include "readFields.H"
#include "tecplotWriter.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoField>
void print(const char* msg, Ostream& os, const PtrList<const GeoField>& flds)
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
    const List<wordRe>& excludePatches  //HashSet<word>& excludePatches
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
        else if (stringOps::match(excludePatches, pp.name()))
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Tecplot binary file format writer"
    );
    timeSelector::addOptions();
    #include "addRegionOption.H"
    argList::addOption
    (
        "fields",
        "names",
        "convert selected fields only. eg, '(p T U)'"
    );
    argList::addOption
    (
        "cellSet",
        "name",
        "restrict conversion to the specified cellSet"
    );
    argList::addOption
    (
        "faceSet",
        "name",
        "restrict conversion to the specified faceSet"
    );
    argList::addBoolOption
    (
        "nearCellValue",
        "output cell value on patches instead of patch value itself"
    );
    argList::addBoolOption
    (
        "noInternal",
        "do not generate file for mesh, only for patches"
    );
    argList::addBoolOption
    (
        "noPointValues",
        "no pointFields"
    );
    argList::addOption
    (
        "excludePatches",
        "patches (wildcards) to exclude"
    );
    argList::addBoolOption
    (
        "noFaceZones",
        "no faceZones"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    const bool doWriteInternal = !args.optionFound("noInternal");
    const bool doFaceZones     = !args.optionFound("noFaceZones");
    const bool nearCellValue = args.optionFound("nearCellValue");
    const bool noPointValues = args.optionFound("noPointValues");

    if (nearCellValue)
    {
        WarningInFunction
            << "Using neighbouring cell value instead of patch value"
            << nl << endl;
    }

    if (noPointValues)
    {
        WarningInFunction
            << "Outputting cell values only" << nl << endl;
    }

    List<wordRe> excludePatches;
    if (args.optionFound("excludePatches"))
    {
        args.optionLookup("excludePatches")() >> excludePatches;

        Info<< "Not including patches " << excludePatches << nl << endl;
    }

    word cellSetName;
    word faceSetName;
    string pltName = runTime.caseName();

    if (args.optionReadIfPresent("cellSet", cellSetName))
    {
        pltName = cellSetName;
    }
    else if (Pstream::parRun())
    {
        // Strip off leading casename, leaving just processor_DDD ending.
        pltName = runTime.caseName();

        string::size_type i = pltName.rfind("processor");

        if (i != string::npos)
        {
            pltName = pltName.substr(i);
        }
    }
    args.optionReadIfPresent("faceSet", faceSetName);

    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createNamedMesh.H"

    // TecplotData/ directory in the case
    fileName fvPath(runTime.path()/"Tecplot360");

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
         || cellSetName.size()
         || faceSetName.size()
         || regionName != polyMesh::defaultRegion
        )
        {
            Info<< "Keeping old tecplot files in " << fvPath << nl << endl;
        }
        else
        {
            Info<< "Deleting old tecplot files in " << fvPath << nl << endl;

            rmDir(fvPath);
        }
    }

    mkDir(fvPath);

    // Mesh wrapper: does subsetting
    meshSubsetHelper meshRef(mesh, meshSubsetHelper::SET, cellSetName);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time: " << runTime.timeName() << endl;

        const word timeDesc = name(timeI); // Foam::name(runTime.timeIndex());

        // Check for new polyMesh/ and update mesh, fvMeshSubset and cell
        // decomposition.
        polyMesh::readUpdateState meshState = meshRef.readUpdate();
        const fvMesh& mesh = meshRef.mesh();

        // TotalNumFaceNodes
        int32_t nFaceNodes = 0;
        forAll(mesh.faces(), facei)
        {
            nFaceNodes += mesh.faces()[facei].size();
        }

        // Read all fields on the new mesh
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Search for list of objects for this time
        IOobjectList objects(mesh, runTime.timeName());

        HashSet<word> selectedFields;
        if (args.optionFound("fields"))
        {
            args.optionLookup("fields")() >> selectedFields;
        }

        // Construct the vol fields (on the original mesh if subsetted)

        PtrList<const volScalarField> vsf;
        readFields(meshRef, meshRef.baseMesh(), objects, selectedFields, vsf);
        print("    volScalarFields            :", Info, vsf);

        PtrList<const volVectorField> vvf;
        readFields(meshRef, meshRef.baseMesh(), objects, selectedFields, vvf);
        print("    volVectorFields            :", Info, vvf);

        PtrList<const volSphericalTensorField> vSphf;
        readFields(meshRef, meshRef.baseMesh(), objects, selectedFields, vSphf);
        print("    volSphericalTensorFields   :", Info, vSphf);

        PtrList<const volSymmTensorField> vSymf;
        readFields(meshRef, meshRef.baseMesh(), objects, selectedFields, vSymf);
        print("    volSymmTensorFields        :", Info, vSymf);

        PtrList<const volTensorField> vtf;
        readFields(meshRef, meshRef.baseMesh(), objects, selectedFields, vtf);
        print("    volTensorFields            :", Info, vtf);


        // Construct pointMesh only if nessecary since constructs edge
        // addressing (expensive on polyhedral meshes)
        if (noPointValues)
        {
            Info<< "    pointScalarFields : switched off"
                << " (\"-noPointValues\" (at your option)\n";
            Info<< "    pointVectorFields : switched off"
                << " (\"-noPointValues\" (at your option)\n";
        }

        PtrList<const pointScalarField> psf;
        PtrList<const pointVectorField> pvf;
        //PtrList<const pointSphericalTensorField> pSphf;
        //PtrList<const pointSymmTensorField> pSymf;
        //PtrList<const pointTensorField> ptf;


        if (!noPointValues)
        {
            //// Add interpolated volFields
            //const volPointInterpolation& pInterp = volPointInterpolation::New
            //(
            //    mesh
            //);
            //
            //label nPsf = psf.size();
            //psf.setSize(nPsf+vsf.size());
            //forAll(vsf, i)
            //{
            //    Info<< "Interpolating " << vsf[i].name() << endl;
            //    tmp<pointScalarField> tvsf(pInterp.interpolate(vsf[i]));
            //    tvsf().rename(vsf[i].name() + "_point");
            //    psf.set(nPsf, tvsf);
            //    nPsf++;
            //}
            //
            //label nPvf = pvf.size();
            //pvf.setSize(vvf.size());
            //forAll(vvf, i)
            //{
            //    Info<< "Interpolating " << vvf[i].name() << endl;
            //    tmp<pointVectorField> tvvf(pInterp.interpolate(vvf[i]));
            //    tvvf().rename(vvf[i].name() + "_point");
            //    pvf.set(nPvf, tvvf);
            //    nPvf++;
            //}

            readFields
            (
                meshRef,
                pointMesh::New(meshRef.baseMesh()),
                objects,
                selectedFields,
                psf
            );
            print("    pointScalarFields          :", Info, psf);

            readFields
            (
                meshRef,
                pointMesh::New(meshRef.baseMesh()),
                objects,
                selectedFields,
                pvf
            );
            print("    pointVectorFields          :", Info, pvf);

            //readFields
            //(
            //    meshRef,
            //    pointMesh::New(meshRef.baseMesh()),
            //    objects,
            //    selectedFields,
            //    pSphf
            //);
            //print("    pointSphericalTensorFields :", Info, pSphf);
            //
            //readFields
            //(
            //    meshRef,
            //    pointMesh::New(meshRef.baseMesh()),
            //    objects,
            //    selectedFields,
            //    pSymf
            //);
            //print("    pointSymmTensorFields      :", Info, pSymf);
            //
            //readFields
            //(
            //    meshRef,
            //    pointMesh::New(meshRef.baseMesh()),
            //    objects,
            //    selectedFields,
            //    ptf
            //);
            //print("    pointTensorFields          :", Info, ptf);
        }
        Info<< endl;


        // Get field names
        // ~~~~~~~~~~~~~~~

        string varNames;
        DynamicList<int32_t> varLocation;

        string cellVarNames;
        DynamicList<int32_t> cellVarLocation;

        // volFields
        tecplotWriter::getTecplotNames
        (
            vsf,
            tecplotWriter::CELL_CENTERED,
            varNames,
            varLocation
        );
        tecplotWriter::getTecplotNames
        (
            vsf,
            tecplotWriter::CELL_CENTERED,
            cellVarNames,
            cellVarLocation
        );

        tecplotWriter::getTecplotNames
        (
            vvf,
            tecplotWriter::CELL_CENTERED,
            varNames,
            varLocation
        );
        tecplotWriter::getTecplotNames
        (
            vvf,
            tecplotWriter::CELL_CENTERED,
            cellVarNames,
            cellVarLocation
        );

        tecplotWriter::getTecplotNames
        (
            vSphf,
            tecplotWriter::CELL_CENTERED,
            varNames,
            varLocation
        );
        tecplotWriter::getTecplotNames
        (
            vSphf,
            tecplotWriter::CELL_CENTERED,
            cellVarNames,
            cellVarLocation
        );

        tecplotWriter::getTecplotNames
        (
            vSymf,
            tecplotWriter::CELL_CENTERED,
            varNames,
            varLocation
        );
        tecplotWriter::getTecplotNames
        (
            vSymf,
            tecplotWriter::CELL_CENTERED,
            cellVarNames,
            cellVarLocation
        );

        tecplotWriter::getTecplotNames
        (
            vtf,
            tecplotWriter::CELL_CENTERED,
            varNames,
            varLocation
        );
        tecplotWriter::getTecplotNames
        (
            vtf,
            tecplotWriter::CELL_CENTERED,
            cellVarNames,
            cellVarLocation
        );


        // pointFields
        tecplotWriter::getTecplotNames
        (
            psf,
            tecplotWriter::NODE_CENTERED,
            varNames,
            varLocation
        );

        tecplotWriter::getTecplotNames
        (
            pvf,
            tecplotWriter::NODE_CENTERED,
            varNames,
            varLocation
        );

        // strandID (= piece id).
        // Gets incremented for every piece of geometry that is output.
        int32_t strandID = 1;

        if (meshState != polyMesh::UNCHANGED)
        {
            if (doWriteInternal)
            {
                // Output mesh and fields
                fileName pltFileName
                (
                    fvPath/pltName
                  + "_"
                  + timeDesc
                  + ".plt"
                );

                const string allVarNames = tecplotWriter::XYZ + " " + varNames;
                DynamicList<int32_t> allVarLocation
                {
                    tecplotWriter::NODE_CENTERED,
                    tecplotWriter::NODE_CENTERED,
                    tecplotWriter::NODE_CENTERED
                };
                allVarLocation.append(varLocation);


                tecplotWriter writer(runTime);
                writer.writeInit
                (
                    runTime.caseName(),
                    allVarNames,
                    pltFileName,
                    tecplotWriter::FILETYPE_FULL
                );

                writer.writePolyhedralZone
                (
                    mesh.name(),        // regionName
                    strandID++,         // strandID
                    mesh,
                    allVarLocation,
                    nFaceNodes
                );

                // Coordinates
                writer.writeField(mesh.points());

                // Write all fields
                writer.writeFields(vsf);
                writer.writeFields(vvf);
                writer.writeFields(vSphf);
                writer.writeFields(vSymf);
                writer.writeFields(vtf);

                writer.writeFields(psf);
                writer.writeFields(pvf);

                writer.writeConnectivity(mesh);
                writer.writeEnd();
            }
        }
        else
        {
            if (doWriteInternal)
            {
                if (timeI == 0)
                {
                    // Output static mesh only
                    fileName pltFileName
                    (
                        fvPath/pltName
                      + "_grid_"
                      + timeDesc
                      + ".plt"
                    );


                    tecplotWriter writer(runTime);
                    writer.writeInit
                    (
                        runTime.caseName(),
                        tecplotWriter::XYZ,
                        pltFileName,
                        tecplotWriter::FILETYPE_GRID
                    );

                    writer.writePolyhedralZone
                    (
                        mesh.name(),        // regionName
                        strandID,           // strandID
                        mesh,
                        List<int32_t>(3, tecplotWriter::NODE_CENTERED),
                        nFaceNodes
                    );

                    // Coordinates
                    writer.writeField(mesh.points());
                    writer.writeConnectivity(mesh);
                    writer.writeEnd();
                }

                // Output solution file
                fileName pltFileName
                (
                    fvPath/pltName
                  + "_"
                  + timeDesc
                  + ".plt"
                );


                tecplotWriter writer(runTime);
                writer.writeInit
                (
                    runTime.caseName(),
                    varNames,
                    pltFileName,
                    tecplotWriter::FILETYPE_SOLUTION
                );

                writer.writePolyhedralZone
                (
                    mesh.name(),        // regionName
                    strandID++,         // strandID
                    mesh,
                    varLocation,
                    0
                );

                // Write all fields
                writer.writeFields(vsf);
                writer.writeFields(vvf);
                writer.writeFields(vSphf);
                writer.writeFields(vSymf);
                writer.writeFields(vtf);

                writer.writeFields(psf);
                writer.writeFields(pvf);

                writer.writeEnd();
            }
        }


        //---------------------------------------------------------------------
        //
        // Write faceSet
        //
        //---------------------------------------------------------------------

        if (faceSetName.size())
        {
            // Load the faceSet
            labelList faceLabels(faceSet(mesh, faceSetName).toc());

            // Filename as if patch with same name.
            mkDir(fvPath/faceSetName);

            fileName patchFileName
            (
                fvPath/faceSetName/faceSetName
              + "_"
              + timeDesc
              + ".plt"
            );

            Info<< "    FaceSet   : " << patchFileName << endl;

            const string allVarNames = tecplotWriter::XYZ + " " + cellVarNames;
            DynamicList<int32_t> allVarLocation
            {
                tecplotWriter::NODE_CENTERED,
                tecplotWriter::NODE_CENTERED,
                tecplotWriter::NODE_CENTERED
            };
            allVarLocation.append(cellVarLocation);


            tecplotWriter writer(runTime);
            writer.writeInit
            (
                runTime.caseName(),
                cellVarNames,
                patchFileName,
                tecplotWriter::FILETYPE_FULL
            );

            const indirectPrimitivePatch ipp
            (
                IndirectList<face>(mesh.faces(), faceLabels),
                mesh.points()
            );

            writer.writePolygonalZone
            (
                faceSetName,
                strandID++,
                ipp,
                allVarLocation
            );

            // Coordinates
            writer.writeField(ipp.localPoints());

            // Write all volfields
            forAll(vsf, i)
            {
                writer.writeField
                (
                    writer.getFaceField
                    (
                        linearInterpolate(vsf[i])(),
                        faceLabels
                    )
                );
            }
            forAll(vvf, i)
            {
                writer.writeField
                (
                    writer.getFaceField
                    (
                        linearInterpolate(vvf[i])(),
                        faceLabels
                    )
                );
            }
            forAll(vSphf, i)
            {
                writer.writeField
                (
                    writer.getFaceField
                    (
                        linearInterpolate(vSphf[i])(),
                        faceLabels
                    )
                );
            }
            forAll(vSymf, i)
            {
                writer.writeField
                (
                    writer.getFaceField
                    (
                        linearInterpolate(vSymf[i])(),
                        faceLabels
                    )
                );
            }
            forAll(vtf, i)
            {
                writer.writeField
                (
                    writer.getFaceField
                    (
                        linearInterpolate(vtf[i])(),
                        faceLabels
                    )
                );
            }
            writer.writeConnectivity(ipp);

            continue;
        }


        //---------------------------------------------------------------------
        //
        // Write patches as multi-zone file
        //
        //---------------------------------------------------------------------

        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        labelList patchIDs(getSelectedPatches(patches, excludePatches));

        mkDir(fvPath/"boundaryMesh");

        fileName patchFileName;

        if (meshRef.useSubMesh())
        {
            patchFileName =
                fvPath/"boundaryMesh"/cellSetName
              + "_"
              + timeDesc
              + ".plt";
        }
        else
        {
            patchFileName =
                fvPath/"boundaryMesh"/"boundaryMesh"
              + "_"
              + timeDesc
              + ".plt";
        }

        Info<< "    Combined patches     : " << patchFileName << endl;

        const string allVarNames = tecplotWriter::XYZ + " " + varNames;
        DynamicList<int32_t> allVarLocation
        {
            tecplotWriter::NODE_CENTERED,
            tecplotWriter::NODE_CENTERED,
            tecplotWriter::NODE_CENTERED
        };
        allVarLocation.append(varLocation);


        tecplotWriter writer(runTime);
        writer.writeInit
        (
            runTime.caseName(),
            allVarNames,
            patchFileName,
            tecplotWriter::FILETYPE_FULL
        );

        forAll(patchIDs, i)
        {
            label patchID = patchIDs[i];
            const polyPatch& pp = patches[patchID];
            // int32_t strandID = 1 + i;

            if (pp.size() > 0)
            {
                Info<< "    Writing patch " << patchID
                    << tab << pp.name()
                    << tab << "strand:" << strandID
                    << nl << endl;

                const indirectPrimitivePatch ipp
                (
                    IndirectList<face>(pp, identity(pp.size())),
                    pp.points()
                );

                writer.writePolygonalZone
                (
                    pp.name(),
                    strandID++,     //strandID,
                    ipp,
                    allVarLocation
                );

                // Coordinates
                writer.writeField(ipp.localPoints());

                // Write all fields
                forAll(vsf, i)
                {
                    writer.writeField
                    (
                        writer.getPatchField
                        (
                            nearCellValue,
                            vsf[i],
                            patchID
                        )
                    );
                }
                forAll(vvf, i)
                {
                    writer.writeField
                    (
                        writer.getPatchField
                        (
                            nearCellValue,
                            vvf[i],
                            patchID
                        )
                    );
                }
                forAll(vSphf, i)
                {
                    writer.writeField
                    (
                        writer.getPatchField
                        (
                            nearCellValue,
                            vSphf[i],
                            patchID
                        )
                    );
                }
                forAll(vSymf, i)
                {
                    writer.writeField
                    (
                        writer.getPatchField
                        (
                            nearCellValue,
                            vSymf[i],
                            patchID
                        )
                    );
                }
                forAll(vtf, i)
                {
                    writer.writeField
                    (
                        writer.getPatchField
                        (
                            nearCellValue,
                            vtf[i],
                            patchID
                        )
                    );
                }

                forAll(psf, i)
                {
                    writer.writeField
                    (
                        psf[i].boundaryField()[patchID].patchInternalField()
                    );
                }
                forAll(pvf, i)
                {
                    writer.writeField
                    (
                        pvf[i].boundaryField()[patchID].patchInternalField()
                    );
                }

                writer.writeConnectivity(ipp);
            }
            else
            {
                Info<< "    Skipping zero sized patch " << patchID
                    << tab << pp.name()
                    << nl << endl;
            }
        }
        writer.writeEnd();
        Info<< endl;



        //---------------------------------------------------------------------
        //
        // Write faceZones as multi-zone file
        //
        //---------------------------------------------------------------------

        const faceZoneMesh& zones = mesh.faceZones();

        if (doFaceZones && !zones.empty())
        {
            mkDir(fvPath/"faceZoneMesh");

            fileName patchFileName;

            if (meshRef.useSubMesh())
            {
                patchFileName =
                    fvPath/"faceZoneMesh"/cellSetName
                  + "_"
                  + timeDesc
                  + ".plt";
            }
            else
            {
                patchFileName =
                    fvPath/"faceZoneMesh"/"faceZoneMesh"
                  + "_"
                  + timeDesc
                  + ".plt";
            }

            Info<< "    FaceZone  : " << patchFileName << endl;

            const string allVarNames = tecplotWriter::XYZ + " " + cellVarNames;
            DynamicList<int32_t> allVarLocation
            {
                tecplotWriter::NODE_CENTERED,
                tecplotWriter::NODE_CENTERED,
                tecplotWriter::NODE_CENTERED
            };
            allVarLocation.append(cellVarLocation);


            tecplotWriter writer(runTime);
            writer.writeInit
            (
                runTime.caseName(),
                allVarNames,
                patchFileName,
                tecplotWriter::FILETYPE_FULL
            );

            forAll(zones, zoneI)
            {
                const faceZone& pp = zones[zoneI];

                if (pp.size() > 0)
                {
                    const indirectPrimitivePatch ipp
                    (
                        IndirectList<face>(mesh.faces(), pp),
                        mesh.points()
                    );

                    writer.writePolygonalZone
                    (
                        pp.name(),
                        strandID++, //1+patchIDs.size()+zoneI,    //strandID,
                        ipp,
                        allVarLocation
                    );

                    // Coordinates
                    writer.writeField(ipp.localPoints());

                    // Write all volfields
                    forAll(vsf, i)
                    {
                        writer.writeField
                        (
                            writer.getFaceField
                            (
                                linearInterpolate(vsf[i])(),
                                pp
                            )
                        );
                    }
                    forAll(vvf, i)
                    {
                        writer.writeField
                        (
                            writer.getFaceField
                            (
                                linearInterpolate(vvf[i])(),
                                pp
                            )
                        );
                    }
                    forAll(vSphf, i)
                    {
                        writer.writeField
                        (
                            writer.getFaceField
                            (
                                linearInterpolate(vSphf[i])(),
                                pp
                            )
                        );
                    }
                    forAll(vSymf, i)
                    {
                        writer.writeField
                        (
                            writer.getFaceField
                            (
                                linearInterpolate(vSymf[i])(),
                                pp
                            )
                        );
                    }
                    forAll(vtf, i)
                    {
                        writer.writeField
                        (
                            writer.getFaceField
                            (
                                linearInterpolate(vtf[i])(),
                                pp
                            )
                        );
                    }

                    writer.writeConnectivity(ipp);
                }
                else
                {
                    Info<< "    Skipping zero sized faceZone " << zoneI
                        << tab << pp.name()
                        << nl << endl;
                }
            }

            writer.writeEnd();
            Info<< endl;
        }


        //---------------------------------------------------------------------
        //
        // Write lagrangian data
        //
        //---------------------------------------------------------------------

        fileNameList cloudDirs
        (
            readDir
            (
                runTime.timePath()/regionPrefix/cloud::prefix,
                fileName::DIRECTORY
            )
        );

        forAll(cloudDirs, cloudI)
        {
            IOobjectList sprayObjs
            (
                mesh,
                runTime.timeName(),
                cloud::prefix/cloudDirs[cloudI]
            );

            IOobject* positionsPtr = sprayObjs.lookup("positions");
            IOobject* coordinatesPtr = sprayObjs.lookup("coordinates");

            if (positionsPtr || coordinatesPtr)
            {
                mkDir(fvPath/cloud::prefix/cloudDirs[cloudI]);

                fileName lagrFileName
                (
                    fvPath/cloud::prefix/cloudDirs[cloudI]/cloudDirs[cloudI]
                  + "_" + timeDesc + ".plt"
                );

                Info<< "    Lagrangian: " << lagrFileName << endl;

                wordList labelNames(sprayObjs.names(labelIOField::typeName));
                Info<< "        labels            :";
                print(Info, labelNames);

                wordList scalarNames(sprayObjs.names(scalarIOField::typeName));
                Info<< "        scalars           :";
                print(Info, scalarNames);

                wordList vectorNames(sprayObjs.names(vectorIOField::typeName));
                Info<< "        vectors           :";
                print(Info, vectorNames);

                //wordList sphereNames
                //(
                //    sprayObjs.names
                //    (
                //        sphericalTensorIOField::typeName
                //    )
                //);
                //Info<< "        spherical tensors :";
                //print(Info, sphereNames);
                //
                //wordList symmNames
                //(
                //    sprayObjs.names
                //    (
                //        symmTensorIOField::typeName
                //    )
                //);
                //Info<< "        symm tensors      :";
                //print(Info, symmNames);
                //
                //wordList tensorNames
                //(
                //    sprayObjs.names(tensorIOField::typeName)
                //);
                //Info<< "        tensors           :";
                //print(Info, tensorNames);


                // Load cloud positions
                passiveParticleCloud parcels(mesh, cloudDirs[cloudI]);

                // Get positions as pointField
                pointField positions(parcels.size());
                label n = 0;
                forAllConstIter(passiveParticleCloud, parcels, elmnt)
                {
                    positions[n++] = elmnt().position();
                }


                string allVarNames = tecplotWriter::XYZ;
                DynamicList<int32_t> allVarLocation
                {
                    tecplotWriter::NODE_CENTERED,
                    tecplotWriter::NODE_CENTERED,
                    tecplotWriter::NODE_CENTERED
                };

                tecplotWriter::getTecplotNames<label>
                (
                    labelNames,
                    tecplotWriter::NODE_CENTERED,
                    allVarNames,
                    allVarLocation
                );

                tecplotWriter::getTecplotNames<scalar>
                (
                    scalarNames,
                    tecplotWriter::NODE_CENTERED,
                    allVarNames,
                    allVarLocation
                );

                tecplotWriter::getTecplotNames<vector>
                (
                    vectorNames,
                    tecplotWriter::NODE_CENTERED,
                    allVarNames,
                    allVarLocation
                );


                tecplotWriter writer(runTime);
                writer.writeInit
                (
                    runTime.caseName(),
                    allVarNames,
                    lagrFileName,
                    tecplotWriter::FILETYPE_FULL
                );

                writer.writeOrderedZone
                (
                    cloudDirs[cloudI],
                    strandID++,     //strandID,
                    parcels.size(),
                    allVarLocation
                );

                // Coordinates
                writer.writeField(positions);

                // labelFields
                forAll(labelNames, i)
                {
                    const word& fieldName = labelNames[i];

                    IOField<label> fld
                    (
                        IOobject
                        (
                            fieldName,
                            mesh.time().timeName(),
                            cloud::prefix/cloudDirs[cloudI],
                            mesh,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE,
                            false
                        )
                    );

                    writer.writeField(fld);
                }

                // scalarFields
                forAll(scalarNames, i)
                {
                    const word& fieldName = scalarNames[i];

                    IOField<scalar> fld
                    (
                        IOobject
                        (
                            fieldName,
                            mesh.time().timeName(),
                            cloud::prefix/cloudDirs[cloudI],
                            mesh,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE,
                            false
                        )
                    );

                    writer.writeField(fld);
                }

                // vectorFields
                forAll(vectorNames, i)
                {
                    const word& fieldName = vectorNames[i];

                    IOField<vector> fld
                    (
                        IOobject
                        (
                            fieldName,
                            mesh.time().timeName(),
                            cloud::prefix/cloudDirs[cloudI],
                            mesh,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE,
                            false
                        )
                    );

                    writer.writeField(fld);
                }

                writer.writeEnd();
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
