/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

InClass
    vtkPV3Foam

\*---------------------------------------------------------------------------*/

#ifndef vtkPV3FoamVolFields_H
#define vtkPV3FoamVolFields_H

// Foam includes
#include "emptyFvPatchField.H"
#include "wallPolyPatch.H"
#include "faceSet.H"
#include "vtkPV3FoamFaceField.H"
#include "vtkPV3FoamPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtkPV3Foam::convertVolFields
(
    const fvMesh& mesh,
    const volPointInterpolation& pInterp,
    const PtrList<PrimitivePatchInterpolation<primitivePatch> >& ppInterpList,
    const IOobjectList& objects,
    vtkMultiBlockDataSet* output
)
{
    // field subset based on type
    IOobjectList fieldObjects
    (
        objects.lookupClass
        (
            GeometricField<Type, fvPatchField, volMesh>::typeName
        )
    );

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    vtkDataArraySelection* regionSelection = reader_->GetRegionSelection();

    forAllIter(IOobjectList, fieldObjects, iter)
    {
        word fieldName = iter()->name();

        GeometricField<Type, fvPatchField, volMesh> tf
        (
            *iter(),
            mesh
        );

        tmp<GeometricField<Type, pointPatchField, pointMesh> > tptf
        (
            pInterp.interpolate(tf)
        );


        //
        // Convert internal mesh - if activated
        //
        for
        (
            int regionId = regionInfoVolume_.start();
            regionId < regionInfoVolume_.end();
            ++regionId
        )
        {
            const label datasetNo = regionDataset_[regionId];

            if (!regionStatus_[regionId] || datasetNo < 0)
            {
                continue;
            }

            convertVolField
            (
                tf,
                output,
                regionInfoVolume_,
                datasetNo,
                superCells_
            );
            convertPointField
            (
                tptf(),
                tf,
                output,
                regionInfoVolume_,
                datasetNo
            );
        }

        //
        // Convert patches - if activated
        //
        for
        (
            int regionId = regionInfoPatches_.start();
            regionId < regionInfoPatches_.end();
            ++regionId
        )
        {
            word patchName = getFirstWord
            (
                regionSelection->GetArrayName(regionId)
            );

            const label datasetNo = regionDataset_[regionId];
            const label patchId = patches.findPatchID(patchName);

            if (!regionStatus_[regionId] || datasetNo < 0 || patchId < 0)
            {
                continue;
            }

            const fvPatchField<Type>& ptf = tf.boundaryField()[patchId];

            if
            (
                isType<emptyFvPatchField<Type> >(ptf)
             ||
                (
                    typeid(patches[patchId]) == typeid(wallPolyPatch)
                 && reader_->GetExtrapolateWalls()
                )
            )
            {
                fvPatch p
                (
                    ptf.patch().patch(),
                    tf.mesh().boundary()
                );

                convertPatchField
                (
                    tf.name(),
                    fvPatchField<Type>(p, tf).patchInternalField()(),
                    output,
                    regionInfoPatches_,
                    datasetNo
                );
            }
            else
            {
                convertPatchField
                (
                    tf.name(),
                    ptf,
                    output,
                    regionInfoPatches_,
                    datasetNo
                );
            }

            convertPatchPointField
            (
                tf.name(),
                ppInterpList[patchId].faceToPointInterpolate(ptf)(),
                output,
                regionInfoPatches_,
                datasetNo
            );
        }

        //
        // Convert cell zones - if activated
        //
        for
        (
            int regionId = regionInfoCellZones_.start();
            regionId < regionInfoCellZones_.end();
            ++regionId
        )
        {
            const label datasetNo = regionDataset_[regionId];

            if (!regionStatus_[regionId] || datasetNo < 0)
            {
                continue;
            }

            if (debug)
            {
                Info<< "wish to convert cellzone: "
                    << getFirstWord(regionSelection->GetArrayName(regionId))
                    << " regionId: " << regionId
                    << " volume field: " << fieldName
                    << endl;
            }

            convertVolField
            (
                tf,
                output,
                regionInfoCellZones_,
                datasetNo,
                zoneSuperCells_[datasetNo]
            );
        }

        //
        // Convert cell sets - if activated
        //
        for
        (
            int regionId = regionInfoCellSets_.start();
            regionId < regionInfoCellSets_.end();
            ++regionId
        )
        {
            const label datasetNo = regionDataset_[regionId];

            if (!regionStatus_[regionId] || datasetNo < 0)
            {
                continue;
            }

            if (debug)
            {
                Info<< "wish to convert cellset: "
                    << getFirstWord(regionSelection->GetArrayName(regionId))
                    << " regionId: " << regionId
                    << " volume field: " << fieldName
                    << endl;
            }

            convertVolField
            (
                tf,
                output,
                regionInfoCellSets_,
                datasetNo,
                csetSuperCells_[datasetNo]
            );
        }

        //
        // Convert face zones - if activated
        //
        for
        (
            int regionId = regionInfoFaceZones_.start();
            regionId < regionInfoFaceZones_.end();
            ++regionId
        )
        {
            word zoneName = getFirstWord
            (
                regionSelection->GetArrayName(regionId)
            );
            const label datasetNo = regionDataset_[regionId];

            if (!regionStatus_[regionId] || datasetNo < 0)
            {
                continue;
            }

            const faceZoneMesh& zMesh = mesh.faceZones();
            const label zoneId = zMesh.findZoneID(zoneName);

            if (zoneId < 0)
            {
                continue;
            }

            if (debug)
            {
                Info<< "wish to convert facezone: "
                    << zoneName
                    << " regionId: " << regionId
                    << " volume field: " << fieldName
                    << endl;
            }

            convertFaceField
            (
                tf,
                output,
                regionInfoFaceZones_,
                datasetNo,
                mesh,
                zMesh[zoneId]
            );
        }

        //
        // Convert face sets - if activated
        //
        for
        (
            int regionId = regionInfoFaceSets_.start();
            regionId < regionInfoFaceSets_.end();
            ++regionId
        )
        {
            const label datasetNo = regionDataset_[regionId];
            if (!regionStatus_[regionId] || datasetNo < 0)
            {
                continue;
            }

            word selectName = getFirstWord
            (
                regionSelection->GetArrayName(regionId)
            );

            if (debug)
            {
                Info<< "wish to convert faceset: " << selectName
                    << " regionId: " << regionId
                    << " volume field: " << fieldName
                    << endl;
            }

            const faceSet fSet(mesh, selectName);

            convertFaceField
            (
                tf,
                output,
                regionInfoFaceSets_,
                datasetNo,
                mesh,
                fSet
            );
        }
    }
}


template<class Type>
void Foam::vtkPV3Foam::convertVolField
(
    const GeometricField<Type, fvPatchField, volMesh>& tf,
    vtkMultiBlockDataSet* output,
    const selectionInfo& selector,
    const label datasetNo,
    labelList& superCells
)
{
    const label nComp = pTraits<Type>::nComponents;

    vtkUnstructuredGrid* vtkmesh = vtkUnstructuredGrid::SafeDownCast
    (
        GetDataSetFromBlock(output, selector, datasetNo)
    );

    vtkFloatArray* celldata = vtkFloatArray::New();
    celldata->SetNumberOfTuples( superCells.size() );
    celldata->SetNumberOfComponents( nComp );
    celldata->Allocate( nComp*superCells.size() );
    celldata->SetName( tf.name().c_str() );

    if (debug)
    {
        Info<< "convert volField: "
            << tf.name()
            << " size = " << tf.size()
            << " nComp=" << nComp
            << " nTuples = " << superCells.size() <<  endl;
    }

    float vec[nComp];
    forAll(superCells, i)
    {
        const Type& t = tf[superCells[i]];
        for (direction d=0; d<nComp; d++)
        {
            vec[d] = component(t, d);
        }

        celldata->InsertTuple(i, vec);
    }

    vtkmesh->GetCellData()->AddArray(celldata);
    celldata->Delete();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
