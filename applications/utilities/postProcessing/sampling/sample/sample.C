/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Sample field data with a choice of interpolation schemes, sampling options
    and write formats.

    interpolationScheme : choice of
        cell            : use cell-centre value only; constant over cells
        cellPoint       : use cell-centre and vertex values
        cellPointFace   : use cell-centre, vertex and face values.

    sample: choice of
        uniform             evenly distributed points on line
        face                one point per face intersection
        midPoint            one point per cell, inbetween two face intersections
        midPointAndFace     combination of face and midPoint

        curve               specified points, not nessecary on line, uses
                            tracking
        cloud               specified points, uses findCell

    writeFormat : choice of
        xmgr
        jplot
        gnuplot
        raw

\*---------------------------------------------------------------------------*/

#include "Pstream.H"
#include "argList.H"
#include "OSspecific.H"

#include "Cloud.H"
#include "passiveParticle.H"
#include "meshSearch.H"
#include "interpolation.H"
#include "volPointInterpolation.H"

#include "writer.H"
#include "sampleSet.H"
#include "volFieldSampler.H"
#include "dictionaryEntry.H"

#include "combineSampleSets.H"
#include "combineSampleValues.H"

using namespace Foam;


template<class Type>
void writeSampleFile
(
    const coordSet& masterSampleSet,
    const PtrList<volFieldSampler<Type> >& masterFields,
    const label setI,
    const fileName& timeDir,
    const word& writeFormat
)
{
    if (masterFields.size() > 0)
    {
        wordList valueSetNames(masterFields.size());
        List<const Field<Type>*> valueSets(masterFields.size());

        forAll(masterFields, fieldI)
        {
            valueSetNames[fieldI] = masterFields[fieldI].name();
            valueSets[fieldI] = &masterFields[fieldI][setI];
        }

        autoPtr<writer<Type> > formatter
        (
            writer<Type>::New(writeFormat)
        );

        fileName fName
        (
            timeDir/formatter().getFileName(masterSampleSet, valueSetNames)
        );
        
        Info<< "Writing fields to " << fName << endl;
        
        formatter().write
        (
            masterSampleSet,
            valueSetNames,
            valueSets,
            OFstream(fName)()
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
#   include "addTimeOptions.H"
#   include "setRootCase.H"

#   include "createTime.H"

    // Get times list
    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
#   include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);

#   include "createMesh.H"

    //
    // Hack: initialize Cloud to initialize the processor table so from
    // now on we can use cloud on single processors only.
    //
    Cloud<passiveParticle> dummyCloud(mesh, IDLList<passiveParticle>());

    // Read control dictionary
    IOdictionary sampleDict
    (
        IOobject
        (
            "sampleDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word interpolationScheme(sampleDict.lookup("interpolationScheme"));
    const wordList fieldNames = sampleDict.lookup("fields");

    //
    // Construct writers
    //

    word writeFormat(sampleDict.lookup("writeFormat"));

    //
    // Construct interpolation dictionary (same interpolation for all fields)
    //

    dictionary interpolationSchemes;

    forAll(fieldNames, fieldI)
    {
        interpolationSchemes.add
        (
            fieldNames[fieldI],
            interpolationScheme
        );
    }

    // Set up interpolation
    autoPtr<pointMesh> pMeshPtr(new pointMesh(mesh));
    autoPtr<volPointInterpolation> pInterpPtr
    (
        new volPointInterpolation(mesh, pMeshPtr())
    );

    // Set up mesh searching
    meshSearch searchEngine(mesh, true);


    fileName samplePath;

    if (Pstream::master())
    {
        if (Pstream::parRun())
        {
            samplePath = runTime.path()/".."/"samples";
        }
        else
        {
            samplePath = runTime.path()/"samples";
        }

        if (exists(samplePath))
        {
            Info<< "Deleting samples/ directory" << endl << endl;
            rmDir(samplePath);
        }
    }

    fileName oldPointsDir("constant");

    for (label i=startTime; i<endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        //
        // Handle geometry/topology changes
        //
        polyMesh::readUpdateState state = mesh.readUpdate();

        if
        (
            state == polyMesh::POINTS_MOVED
         || state == polyMesh::TOPO_CHANGE
        )
        {
            // Geometry and topology changes            
            searchEngine.correct();

            pMeshPtr.reset(new pointMesh(mesh));

            pInterpPtr.reset(new volPointInterpolation(mesh, pMeshPtr()));
        }


        //
        // Construct sampling point generators
        //

        PtrList<sampleSet> sampleSets
        (
            sampleDict.lookup("sampleSets"),
            sampleSet::iNew(mesh, searchEngine)
        );
        if (sampleSets.size() < 1)
        {
            FatalErrorIn(args.executable())
                << "No sampleSets provided in sampleDict"
                << exit(FatalError);
        }

        // Storage for interpolated values
        PtrList<volFieldSampler<scalar> > sampledScalarFields
        (
            fieldNames.size()
        );
        PtrList<volFieldSampler<vector> > sampledVectorFields
        (
            fieldNames.size()
        );
        PtrList<volFieldSampler<sphericalTensor> > sampledSphericalTensorFields
        (
            fieldNames.size()
        );
        PtrList<volFieldSampler<symmTensor> > sampledSymmTensorFields
        (
            fieldNames.size()
        );
        PtrList<volFieldSampler<tensor> > sampledTensorFields
        (
            fieldNames.size()
        );

        //
        // Do actual interpolation
        //

        label nScalarFields = 0;
        label nVectorFields = 0;
        label nSphericalTensorFields = 0;
        label nSymmTensorFields = 0;
        label nTensorFields = 0;

        forAll(fieldNames, fieldI)
        {
            const word& fieldName = fieldNames[fieldI];

            IOobject fieldHeader
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            );

            // Determine number of processor actually having this field
            label fieldFound = (fieldHeader.headerOk() ? 1 : 0);
            reduce(fieldFound, sumOp<label>());

            if (fieldFound == Pstream::nProcs())
            {
                if 
                (
                    fieldHeader.headerClassName() == volScalarField::typeName
                )
                {
                    Info<< "Sampling " << fieldHeader.headerClassName()
                        << ' ' << fieldName << endl;

                    volScalarField sField(fieldHeader, mesh);

                    sampledScalarFields.set
                    (
                        nScalarFields,
                        new volFieldSampler<scalar>
                        (
                            pInterpPtr(),
                            interpolationSchemes,
                            sField,
                            sampleSets
                        )
                    );

                    nScalarFields++;
                }
                else if
                (
                    fieldHeader.headerClassName() == volVectorField::typeName
                )
                {
                    Info<< "Sampling " << fieldHeader.headerClassName()
                        << ' ' << fieldName << endl;

                    volVectorField vField(fieldHeader, mesh);

                    sampledVectorFields.set
                    (
                        nVectorFields,
                        new volFieldSampler<vector>
                        (
                            pInterpPtr(),
                            interpolationSchemes,
                            vField,
                            sampleSets
                        )
                    );

                    nVectorFields++;
                }
                else if
                (
                    fieldHeader.headerClassName()
                 == volSphericalTensorField::typeName
                )
                {
                    Info<< "Sampling " << fieldHeader.headerClassName()
                        << ' ' << fieldName << endl;

                    volSphericalTensorField tField(fieldHeader, mesh);

                    sampledSphericalTensorFields.set
                    (
                        nSphericalTensorFields,
                        new volFieldSampler<sphericalTensor>
                        (
                            pInterpPtr(),
                            interpolationSchemes,
                            tField,
                            sampleSets
                        )
                    );

                    nSphericalTensorFields++;
                }
                else if
                (
                    fieldHeader.headerClassName()
                 == volSymmTensorField::typeName
                )
                {
                    Info<< "Sampling " << fieldHeader.headerClassName()
                        << ' ' << fieldName << endl;

                    volSymmTensorField tField(fieldHeader, mesh);

                    sampledSymmTensorFields.set
                    (
                        nSymmTensorFields,
                        new volFieldSampler<symmTensor>
                        (
                            pInterpPtr(),
                            interpolationSchemes,
                            tField,
                            sampleSets
                        )
                    );

                    nSymmTensorFields++;
                }
                else if
                (
                    fieldHeader.headerClassName() == volTensorField::typeName
                )
                {
                    Info<< "Sampling " << fieldHeader.headerClassName()
                        << ' ' << fieldName << endl;

                    volTensorField tField(fieldHeader, mesh);

                    sampledTensorFields.set
                    (
                        nTensorFields,
                        new volFieldSampler<tensor>
                        (
                            pInterpPtr(),
                            interpolationSchemes,
                            tField,
                            sampleSets
                        )
                    );

                    nTensorFields++;
                }
            }
            else if (fieldFound != 0)
            {
                FatalErrorIn(args.executable())
                    << "Did not find field " << fieldName
                    << " on all processors" << exit(FatalError);
            }
            else if (fieldName.find('(') != string::npos)
            {
                if (fieldName.find("component") != string::npos)
                {
                    string baseFieldName(fieldName(0, fieldName.find('.')));

                    IOobject fieldHeader
                    (
                        baseFieldName,
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ
                    );

                    // Determine number of processor actually having this field
                    label fieldFound = (fieldHeader.headerOk() ? 1 : 0);
                    reduce(fieldFound, sumOp<label>());

                    if (fieldFound == Pstream::nProcs())
                    {
                        if
                        (
                            fieldHeader.headerClassName()
                         == volVectorField::typeName
                        )
                        {
                            size_t cmptPos(fieldName.find_last_of("012"));

                            if (cmptPos == string::npos)
                            {
                                FatalErrorIn(args.executable())
                                    << "cannot find component index for "
                                    << fieldName << exit(FatalError);
                            }

                            direction cmpt =
                                atoi(string(fieldName[cmptPos]).c_str());

                            Info<< "Sampling " << fieldHeader.headerClassName()
                                << ' ' << fieldName << endl;

                            volVectorField vField(fieldHeader, mesh);

                            sampledScalarFields.set
                            (
                                nScalarFields,
                                new volFieldSampler<scalar>
                                (
                                    pInterpPtr(),
                                    interpolationSchemes,
                                    vField.component(cmpt)(),
                                    sampleSets
                                )
                            );

                            nScalarFields++;
                        }
                        else if
                        (
                            fieldHeader.headerClassName()
                         == volSymmTensorField::typeName
                        )
                        {
                            size_t cmptPos(fieldName.find_last_of("0123456"));

                            if (cmptPos == string::npos)
                            {
                                FatalErrorIn(args.executable())
                                    << "cannot find component index for "
                                    << fieldName
                                    << exit(FatalError);
                            }

                            direction cmpt =
                                atoi(string(fieldName[cmptPos]).c_str());

                            Info<< "Sampling " << fieldHeader.headerClassName()
                                << ' ' << fieldName << endl;

                            volSymmTensorField tField(fieldHeader, mesh);

                            sampledScalarFields.set
                            (
                                nScalarFields,
                                new volFieldSampler<scalar>
                                (
                                    pInterpPtr(),
                                    interpolationSchemes,
                                    tField.component(cmpt)(),
                                    sampleSets
                                )
                            );

                            nScalarFields++;
                        }
                        else if
                        (
                            fieldHeader.headerClassName()
                         == volTensorField::typeName
                        )
                        {
                            size_t cmptPos(fieldName.find_last_of("012345678"));

                            if (cmptPos == string::npos)
                            {
                                FatalErrorIn(args.executable())
                                    << "cannot find component index for "
                                    << fieldName
                                    << exit(FatalError);
                            }

                            direction cmpt =
                                atoi(string(fieldName[cmptPos]).c_str());

                            Info<< "Sampling " << fieldHeader.headerClassName()
                                << ' ' << fieldName << endl;

                            volTensorField tField(fieldHeader, mesh);

                            sampledScalarFields.set
                            (
                                nScalarFields,
                                new volFieldSampler<scalar>
                                (
                                    pInterpPtr(),
                                    interpolationSchemes,
                                    tField.component(cmpt)(),
                                    sampleSets
                                )
                            );

                            nScalarFields++;
                        }
                        else
                        {
                            FatalErrorIn(args.executable())
                                << "component function not supported for field "
                                << fieldName << " of type "
                                << fieldHeader.headerClassName()
                                << exit(FatalError);
                        }
                    }
                    else if (fieldFound != 0)
                    {
                        FatalErrorIn(args.executable())
                            << "Did not find field " << baseFieldName
                            << " on all processors" << exit(FatalError);
                    }

                }
                else if (fieldName.find("mag") != string::npos)
                {
                    string baseFieldName
                    (
                        fieldName(fieldName.find('(') + 1,
                        fieldName.find(')') - fieldName.find('(') - 1)
                    );

                    IOobject fieldHeader
                    (
                        baseFieldName,
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ
                    );

                    // Determine number of processor actually having this field
                    label fieldFound = (fieldHeader.headerOk() ? 1 : 0);
                    reduce(fieldFound, sumOp<label>());


                    if (fieldFound == Pstream::nProcs())
                    {
                        if
                        (
                            fieldHeader.headerClassName()
                         == volScalarField::typeName
                        )
                        {
                            Info<< "Sampling " << fieldHeader.headerClassName()
                                << ' ' << fieldName << endl;

                            volScalarField sField(fieldHeader, mesh);

                            sampledScalarFields.set
                            (
                                nScalarFields,
                                new volFieldSampler<scalar>
                                (
                                    pInterpPtr(),
                                    interpolationSchemes,
                                    mag(sField),
                                    sampleSets
                                )
                            );

                            nScalarFields++;
                        }
                        else if
                        (
                            fieldHeader.headerClassName()
                         == volVectorField::typeName
                        )
                        {
                            Info<< "Sampling " << fieldHeader.headerClassName()
                                << ' ' << fieldName << endl;

                            volVectorField vField(fieldHeader, mesh);

                            sampledScalarFields.set
                            (
                                nScalarFields,
                                new volFieldSampler<scalar>
                                (
                                    pInterpPtr(),
                                    interpolationSchemes,
                                    mag(vField),
                                    sampleSets
                                )
                            );

                            nScalarFields++;
                        }
                        else if
                        (
                            fieldHeader.headerClassName()
                         == volSphericalTensorField::typeName
                        )
                        {
                            Info<< "Sampling " << fieldHeader.headerClassName()
                                << ' ' << fieldName << endl;

                            volSphericalTensorField tField(fieldHeader, mesh);

                            sampledScalarFields.set
                            (
                                nScalarFields,
                                new volFieldSampler<scalar>
                                (
                                    pInterpPtr(),
                                    interpolationSchemes,
                                    mag(tField),
                                    sampleSets
                                )
                            );

                            nScalarFields++;
                        }
                        else if
                        (
                            fieldHeader.headerClassName()
                         == volSymmTensorField::typeName
                        )
                        {
                            Info<< "Sampling " << fieldHeader.headerClassName()
                                << ' ' << fieldName << endl;

                            volSymmTensorField tField(fieldHeader, mesh);

                            sampledScalarFields.set
                            (
                                nScalarFields,
                                new volFieldSampler<scalar>
                                (
                                    pInterpPtr(),
                                    interpolationSchemes,
                                    mag(tField),
                                    sampleSets
                                )
                            );

                            nScalarFields++;
                        }
                        else if
                        (
                            fieldHeader.headerClassName()
                         == volTensorField::typeName
                        )
                        {
                            Info<< "Sampling " << fieldHeader.headerClassName()
                                << ' ' << fieldName << endl;

                            volTensorField tField(fieldHeader, mesh);

                            sampledScalarFields.set
                            (
                                nScalarFields,
                                new volFieldSampler<scalar>
                                (
                                    pInterpPtr(),
                                    interpolationSchemes,
                                    mag(tField),
                                    sampleSets
                                )
                            );

                            nScalarFields++;
                        }
                        else
                        {
                            FatalErrorIn(args.executable())
                                << "mag function not supported for field "
                                << fieldName << " of type "
                                << fieldHeader.headerClassName()
                                << exit(FatalError);
                        }
                    }
                    else if (fieldFound != 0)
                    {
                        FatalErrorIn(args.executable())
                            << "Did not find field " << baseFieldName
                            << " on all processors" << exit(FatalError);
                    }
                }
            }
        }

        // Set the sampledFields to the correct size
        sampledScalarFields.setSize(nScalarFields);
        sampledVectorFields.setSize(nVectorFields);
        sampledSphericalTensorFields.setSize(nSphericalTensorFields);
        sampledSymmTensorFields.setSize(nSymmTensorFields);
        sampledTensorFields.setSize(nTensorFields);

        //
        // Now we have all results
        //   - sampleSets       : list of all sampling sets
        //   - sampledXXXFields : list of all sampled fields
        //

        // Combine sampleSets from processors. Sort by curveDist. Return
        // ordering in indexSets.
        // Note: only master results are valid

        PtrList<coordSet> masterSampleSets(sampleSets.size());
        labelListList indexSets(sampleSets.size());
        combineSampleSets(sampleSets, masterSampleSets, indexSets);


        // Combine sampled fields from processors.
        // Note: only master results are valid

        PtrList<volFieldSampler<scalar> > masterScalarFields
        (
            sampledScalarFields.size()
        );
        combineSampleValues(sampledScalarFields, indexSets, masterScalarFields);

        PtrList<volFieldSampler<vector> > masterVectorFields
        (
            sampledVectorFields.size()
        );
        combineSampleValues(sampledVectorFields, indexSets, masterVectorFields);

        PtrList<volFieldSampler<sphericalTensor> > masterSphericalTensorFields
        (
            sampledSphericalTensorFields.size()
        );
        combineSampleValues
        (
            sampledSphericalTensorFields,
            indexSets,
            masterSphericalTensorFields
        );

        PtrList<volFieldSampler<symmTensor> > masterSymmTensorFields
        (
            sampledSymmTensorFields.size()
        );
        combineSampleValues
        (
            sampledSymmTensorFields,
            indexSets,
            masterSymmTensorFields
        );

        PtrList<volFieldSampler<tensor> > masterTensorFields
        (
            sampledTensorFields.size()
        );
        combineSampleValues(sampledTensorFields, indexSets, masterTensorFields);


        //
        // Write each set, each Field type (scalar/vector/tensor) to separate
        // file.
        //

        if (Pstream::master())
        {
            fileName timeDir(samplePath/runTime.timeName());

            // Mirror the time structure under "samples"
            mkDir(timeDir);

            forAll(masterSampleSets, setI)
            {
                writeSampleFile
                (
                    masterSampleSets[setI],
                    masterScalarFields,
                    setI,
                    timeDir,
                    writeFormat
                );

                writeSampleFile
                (
                    masterSampleSets[setI],
                    masterVectorFields,
                    setI,
                    timeDir,
                    writeFormat
                );

                writeSampleFile
                (
                    masterSampleSets[setI],
                    masterSphericalTensorFields,
                    setI,
                    timeDir,
                    writeFormat
                );

                writeSampleFile
                (
                    masterSampleSets[setI],
                    masterSymmTensorFields,
                    setI,
                    timeDir,
                    writeFormat
                );

                writeSampleFile
                (
                    masterSampleSets[setI],
                    masterTensorFields,
                    setI,
                    timeDir,
                    writeFormat
                );
            }

            Info<< endl;
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
