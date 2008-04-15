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

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "OSspecific.H"
#include "meshSearch.H"
#include "interpolation.H"
#include "volPointInterpolation.H"
#include "cuttingPlane.H"
#include "OFstream.H"
#include "Pstream.H"
#include "ListListOps.H"
#include "Cloud.H"
#include "passiveParticle.H"
#include "mergePoints.H"

#include "fieldsCache.H"
#include "surface.H"
#include "surfaceWriter.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Used to offset faces in Pstream::combineOffset
namespace Foam
{

template <>
class offsetOp<face>
{

public:

    face operator()
    (
        const face& x,
        const label offset
    ) const
    {
        face result(x.size());

        forAll(x, xI)
        {
            result[xI] = x[xI] + offset;
        }
        return result;
    }
};

}


void mergePoints
(
    const polyMesh& mesh,
    const scalar mergeTol,
    List<pointField>& allPoints,
    List<faceList>& allFaces,
    labelListList& allOldToNew
)
{
    const boundBox& bb = mesh.globalData().bb();

    scalar mergeDim = mergeTol * mag(bb.max() - bb.min());

    Info<< nl << "Merging all points within " << mergeDim << " meter." << endl;

    allOldToNew.setSize(allPoints.size());

    forAll(allPoints, surfaceI)
    {
        pointField newPoints;
        labelList oldToNew;

        bool hasMerged = mergePoints
        (
            allPoints[surfaceI],
            mergeDim,
            false,                  // verbosity
            oldToNew,
            newPoints
        );

        if (hasMerged)
        {
            // Copy points
            allPoints[surfaceI].transfer(newPoints);

            // Store point mapping
            allOldToNew[surfaceI].transfer(oldToNew);

            // Relabel faces.
            faceList& faces = allFaces[surfaceI];

            forAll(faces, faceI)
            {
                inplaceRenumber(allOldToNew[surfaceI], faces[faceI]);
            }

            Info<< "For surface " << surfaceI << " merged from "
                << allOldToNew[surfaceI].size() << " points down to "
                << allPoints[surfaceI].size() << " points." << endl;
        }
    }
}


template<class T>
void renumberData
(
    const labelList& oldToNewPoints,
    const label newSize,
    Field<T>& values
)
{
    if (oldToNewPoints.size() == values.size())
    {
        inplaceReorder(oldToNewPoints, values);
        values.setSize(newSize);
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
    // now on we can use cloud (in meshSearch) on single processors only.
    //
    Cloud<passiveParticle> dummyCloud(mesh, IDLList<passiveParticle>());

    //
    // Read control dictionary
    //

    IOdictionary sampleDict
    (
        IOobject
        (
            "sampleSurfaceDict",
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

    word writeFormat(sampleDict.lookup("surfaceFormat"));

    autoPtr<surfaceWriter<scalar> > scalarFormatter
    (
        surfaceWriter<scalar>::New(writeFormat)
    );
    autoPtr<surfaceWriter<vector> > vectorFormatter
    (
        surfaceWriter<vector>::New(writeFormat)
    );
    autoPtr<surfaceWriter<sphericalTensor> > sphericalTensorFormatter
    (
        surfaceWriter<sphericalTensor>::New(writeFormat)
    );
    autoPtr<surfaceWriter<symmTensor> > symmTensorFormatter
    (
        surfaceWriter<symmTensor>::New(writeFormat)
    );
    autoPtr<surfaceWriter<tensor> > tensorFormatter
    (
        surfaceWriter<tensor>::New(writeFormat)
    );

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

    fileName samplePath;

    if (Pstream::parRun())
    {
        samplePath = runTime.path()/".."/"sampleSurfaces";
    }
    else
    {
        samplePath = runTime.path()/"sampleSurfaces";
    }

    if (Pstream::master() && exists(samplePath))
    {
        Info<< "Deleting sampleSurfaces/ directory" << endl << endl;
        rmDir(samplePath);
    }

    // Set up interpolation
    autoPtr<pointMesh> pMeshPtr(new pointMesh(mesh));
    autoPtr<volPointInterpolation> pInterpPtr
    (
        new volPointInterpolation(mesh, pMeshPtr())
    );

    // Set up mesh searching
    meshSearch searchEngine(mesh, true);

    // Create sample surfaces
    PtrList<surface> surfaces
    (
        sampleDict.lookup("surfaces"),
        surface::iNew(mesh, searchEngine)
    );
    Info<< endl;


    fileName oldPointsDir("constant");


    for (label i=startTime; i<endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        //
        // Handle geometry/topology changes
        //
        polyMesh::readUpdateState state = mesh.readUpdate();

        bool meshChanged = false;

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

            meshChanged = true;
        }

        if (Pstream::master())
        {
            Info<< "Creating directory " << samplePath/runTime.timeName()
                << endl << endl;

            mkDir(samplePath/runTime.timeName());
        }


        //
        // Cache used fields and interpolators
        //

        fieldsCache<scalar> scalarCache;
        fieldsCache<vector> vectorCache;
        fieldsCache<sphericalTensor> sphericalTensorCache;
        fieldsCache<symmTensor> symmTensorCache;
        fieldsCache<tensor> tensorCache;

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

            if (fieldHeader.headerOk())
            {
                if 
                (
                    fieldHeader.headerClassName() == volScalarField::typeName
                )
                {
                    Info<< "Loading " << fieldHeader.headerClassName()
                        << ' ' << fieldName << endl;

                    volScalarField* fieldPtr
                    (
                        new volScalarField
                        (
                            fieldHeader,
                            mesh
                        )
                    );

                    scalarCache.insert(fieldName, fieldPtr);
                }

                else if 
                (
                    fieldHeader.headerClassName() == volVectorField::typeName
                )
                {
                    Info<< "Loading " << fieldHeader.headerClassName()
                        << ' ' << fieldName << endl;

                    volVectorField* fieldPtr
                    (
                        new volVectorField
                        (
                            fieldHeader,
                            mesh
                        )
                    );

                    vectorCache.insert(fieldName, fieldPtr);
                }

                else if 
                (
                    fieldHeader.headerClassName()
                 == volSphericalTensorField::typeName
                )
                {
                    Info<< "Loading " << fieldHeader.headerClassName()
                        << ' ' << fieldName << endl;

                    volSphericalTensorField* fieldPtr
                    (
                        new volSphericalTensorField
                        (
                            fieldHeader,
                            mesh
                        )
                    );

                    sphericalTensorCache.insert(fieldName, fieldPtr);
                }

                else if 
                (
                    fieldHeader.headerClassName()
                 == volSymmTensorField::typeName
                )
                {
                    Info<< "Loading " << fieldHeader.headerClassName()
                        << ' ' << fieldName << endl;

                    volSymmTensorField* fieldPtr
                    (
                        new volSymmTensorField
                        (
                            fieldHeader,
                            mesh
                        )
                    );

                    symmTensorCache.insert(fieldName, fieldPtr);
                }

                else if 
                (
                    fieldHeader.headerClassName() == volTensorField::typeName
                )
                {
                    Info<< "Loading " << fieldHeader.headerClassName()
                        << ' ' << fieldName << endl;

                    volTensorField* fieldPtr
                    (
                        new volTensorField
                        (
                            fieldHeader,
                            mesh
                        )
                    );

                    tensorCache.insert(fieldName, fieldPtr);
                }
            }
        }



        //
        // Now we have fields cached and know whether mesh has changed.
        // Update surfaces with this information.
        //

        forAll(surfaces, surfaceI)
        {
            surfaces[surfaceI].correct
            (
                meshChanged,
                pInterpPtr,
                interpolationSchemes,
                scalarCache
            );
        }

        
        //
        // Combine surfaces onto master (bug:should not be redone for static
        // meshes)
        //

        List<pointField> allPoints(surfaces.size());
        List<faceList> allFaces(surfaces.size());

        forAll(surfaces, surfaceI)
        {
            // Collect points from all processors
            List<pointField> gatheredPoints(Pstream::nProcs());
            gatheredPoints[Pstream::myProcNo()] = surfaces[surfaceI].points();
            Pstream::gatherList(gatheredPoints);

            if (Pstream::master())
            {
                allPoints[surfaceI] = ListListOps::combine<pointField>
                (
                    gatheredPoints,
                    accessOp<pointField>()
                );
            }

            // Collect faces from all processors and renumber using sizes of
            // gathered points
            List<faceList> gatheredFaces(Pstream::nProcs());
            gatheredFaces[Pstream::myProcNo()] = surfaces[surfaceI].faces();
            Pstream::gatherList(gatheredFaces);

            if (Pstream::master())
            {
                allFaces[surfaceI] = static_cast<const faceList&>
                (
                    ListListOps::combineOffset<faceList>
                    (
                        gatheredFaces,
                        ListListOps::subSizes
                        (
                            gatheredPoints,
                            accessOp<pointField>()
                        ),
                        accessOp<faceList>(),
                        offsetOp<face>()
                    )
                );
            }
        }

        // Merge close points (1E-10 of mesh bounding box)
        labelListList allOldToNewPoints;
        mergePoints(mesh, 1E-10, allPoints, allFaces, allOldToNewPoints);


        //
        // Do actual interpolation
        //

        forAll(fieldNames, fieldI)
        {
            Info<< endl;

            const word& fieldName = fieldNames[fieldI];

            // Scalar fields

            if (scalarCache.found(fieldName))
            {
                Info<< "Sampling volScalarField " << fieldName << endl << endl;

                forAll(surfaces, surfaceI)
                {
                    const surface& s = surfaces[surfaceI];

                    scalarField values
                    (
                        s.interpolate
                        (
                            fieldName,
                            scalarCache,
                            pInterpPtr,
                            interpolationSchemes
                        )
                    );

                    // Collect values from all processors
                    List<scalarField> gatheredValues(Pstream::nProcs());
                    gatheredValues[Pstream::myProcNo()] = values;
                    Pstream::gatherList(gatheredValues);

                    if (Pstream::master())
                    {
                        // Combine values into single field
                        scalarField allValues
                        (
                            ListListOps::combine<scalarField>
                            (
                                gatheredValues,
                                accessOp<scalarField>()
                            )
                        );

                        // Renumber (point data) to correspond to merged points
                        renumberData
                        (
                            allOldToNewPoints[surfaceI],
                            allPoints[surfaceI].size(),
                            allValues
                        );

                        // Write to time directory under sampleSurfaces/
                        scalarFormatter().write
                        (
                            samplePath,
                            runTime.timeName(),
                            s.name(),
                            allPoints[surfaceI],
                            allFaces[surfaceI],
                            fieldName,
                            allValues
                        );
                    }
                }
            }

            // Vector fields

            else if (vectorCache.found(fieldName))
            {
                Info<< "Sampling volVectorField " << fieldName << endl << endl;

                forAll(surfaces, surfaceI)
                {
                    const surface& s = surfaces[surfaceI];

                    vectorField values
                    (
                        s.interpolate
                        (
                            fieldName,
                            vectorCache,
                            pInterpPtr,
                            interpolationSchemes
                        )
                    );

                    // Collect values from all processors
                    List<vectorField> gatheredValues(Pstream::nProcs());
                    gatheredValues[Pstream::myProcNo()] = values;
                    Pstream::gatherList(gatheredValues);

                    if (Pstream::master())
                    {
                        vectorField allValues
                        (
                            ListListOps::combine<vectorField>
                            (
                                gatheredValues,
                                accessOp<vectorField>()
                            )
                        );

                        // Renumber (point data) to correspond to merged points
                        renumberData
                        (
                            allOldToNewPoints[surfaceI],
                            allPoints[surfaceI].size(),
                            allValues
                        );

                        // Write to time directory under sampleSurfaces/
                        vectorFormatter().write
                        (
                            samplePath,
                            runTime.timeName(),
                            s.name(),
                            allPoints[surfaceI],
                            allFaces[surfaceI],
                            fieldName,
                            allValues
                        );
                    }
                }
            }

            // SphericalTensor fields

            else if (sphericalTensorCache.found(fieldName))
            {
                Info<< "Sampling volSphericalTensorField "
                    << fieldName << endl << endl;

                forAll(surfaces, surfaceI)
                {
                    const surface& s = surfaces[surfaceI];

                    sphericalTensorField values
                    (
                        s.interpolate
                        (
                            fieldName,
                            sphericalTensorCache,
                            pInterpPtr,
                            interpolationSchemes
                        )
                    );

                    // Collect values from all processors
                    List<sphericalTensorField>
                        gatheredValues(Pstream::nProcs());

                    gatheredValues[Pstream::myProcNo()] = values;
                    Pstream::gatherList(gatheredValues);

                    if (Pstream::master())
                    {
                        sphericalTensorField allValues
                        (
                            ListListOps::combine<sphericalTensorField>
                            (
                                gatheredValues,
                                accessOp<sphericalTensorField>()
                            )
                        );

                        // Renumber (point data) to correspond to merged points
                        renumberData
                        (
                            allOldToNewPoints[surfaceI],
                            allPoints[surfaceI].size(),
                            allValues
                        );

                        // Write to time directory under sampleSurfaces/
                        sphericalTensorFormatter().write
                        (
                            samplePath,
                            runTime.timeName(),
                            s.name(),
                            allPoints[surfaceI],
                            allFaces[surfaceI],
                            fieldName,
                            allValues
                        );
                    }
                }
            }

            // SymmTensor fields

            else if (symmTensorCache.found(fieldName))
            {
                Info<< "Sampling volSymmTensorField "
                    << fieldName << endl << endl;

                forAll(surfaces, surfaceI)
                {
                    const surface& s = surfaces[surfaceI];

                    symmTensorField values
                    (
                        s.interpolate
                        (
                            fieldName,
                            symmTensorCache,
                            pInterpPtr,
                            interpolationSchemes
                        )
                    );

                    // Collect values from all processors
                    List<symmTensorField> gatheredValues(Pstream::nProcs());
                    gatheredValues[Pstream::myProcNo()] = values;
                    Pstream::gatherList(gatheredValues);

                    if (Pstream::master())
                    {
                        symmTensorField allValues
                        (
                            ListListOps::combine<symmTensorField>
                            (
                                gatheredValues,
                                accessOp<symmTensorField>()
                            )
                        );

                        // Renumber (point data) to correspond to merged points
                        renumberData
                        (
                            allOldToNewPoints[surfaceI],
                            allPoints[surfaceI].size(),
                            allValues
                        );

                        // Write to time directory under sampleSurfaces/
                        symmTensorFormatter().write
                        (
                            samplePath,
                            runTime.timeName(),
                            s.name(),
                            allPoints[surfaceI],
                            allFaces[surfaceI],
                            fieldName,
                            allValues
                        );
                    }
                }
            }

            // Tensor fields

            else if (tensorCache.found(fieldName))
            {
                Info<< "Sampling volTensorField " << fieldName << endl << endl;

                forAll(surfaces, surfaceI)
                {
                    const surface& s = surfaces[surfaceI];

                    tensorField values
                    (
                        s.interpolate
                        (
                            fieldName,
                            tensorCache,
                            pInterpPtr,
                            interpolationSchemes
                        )
                    );

                    // Collect values from all processors
                    List<tensorField> gatheredValues(Pstream::nProcs());
                    gatheredValues[Pstream::myProcNo()] = values;
                    Pstream::gatherList(gatheredValues);

                    if (Pstream::master())
                    {
                        tensorField allValues
                        (
                            ListListOps::combine<tensorField>
                            (
                                gatheredValues,
                                accessOp<tensorField>()
                            )
                        );

                        // Renumber (point data) to correspond to merged points
                        renumberData
                        (
                            allOldToNewPoints[surfaceI],
                            allPoints[surfaceI].size(),
                            allValues
                        );

                        // Write to time directory under sampleSurfaces/
                        tensorFormatter().write
                        (
                            samplePath,
                            runTime.timeName(),
                            s.name(),
                            allPoints[surfaceI],
                            allFaces[surfaceI],
                            fieldName,
                            allValues
                        );
                    }
                }
            }
        }

        Info<< endl;
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
