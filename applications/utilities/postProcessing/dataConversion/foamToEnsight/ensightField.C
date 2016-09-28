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

\*---------------------------------------------------------------------------*/

#include "ensightFile.H"
#include "ensightField.H"
#include "fvMesh.H"
#include "volFields.H"
#include "OFstream.H"
#include "IOmanip.H"
#include "volPointInterpolation.H"
#include "ensightBinaryStream.H"
#include "ensightAsciiStream.H"
#include "globalIndex.H"
#include "ensightPTraits.H"
#include "zeroGradientFvPatchField.H"

using namespace Foam;

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
volField
(
    const fvMeshSubset& meshSubsetter,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    if (meshSubsetter.hasSubMesh())
    {
        tmp<GeometricField<Type, fvPatchField, volMesh>> tfld
        (
            meshSubsetter.interpolate(vf)
        );
        tfld.ref().checkOut();
        tfld.ref().rename(vf.name());
        return tfld;
    }
    else
    {
        return vf;
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
volField
(
    const fvMeshSubset& meshSubsetter,
    const typename GeometricField
    <
        Type,
        fvPatchField,
        volMesh
    >::DimensionedInternalField& df
)
{
    // Construct volField (with zeroGradient) from dimensioned field

    IOobject io(df);
    io.readOpt()  = IOobject::NO_READ;
    io.writeOpt() = IOobject::NO_WRITE;
    io.registerObject() = false;

    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            io,
            df.mesh(),
            dimensioned<Type>("0", df.dimensions(), Zero),
            zeroGradientFvPatchField<Type>::typeName
        )
    );
    tvf.ref().internalField() = df;
    tvf.ref().correctBoundaryConditions();

    if (meshSubsetter.hasSubMesh())
    {
        const GeometricField<Type, fvPatchField, volMesh>& vf = tvf();

        tmp<GeometricField<Type, fvPatchField, volMesh>> tfld
        (
            meshSubsetter.interpolate(vf)
        );
        tfld.ref().checkOut();
        tfld.ref().rename(vf.name());
        return tfld;
    }
    else
    {
        return tvf;
    }
}


//template<class Container>
//void readAndConvertField
//(
//    const fvMeshSubset& meshSubsetter,
//    const IOobject& io,
//    const fvMesh& mesh,
//    const ensightMesh& eMesh,
//    const fileName& dataDir,
//    const label timeIndex,
//    const bool nodeValues,
//    Ostream& ensightCaseFile
//)
//{
//    Container fld(io, mesh);
//    ensightField<typename Container::value_type>
//    (
//        volField<typename Container::value_type>(meshSubsetter, fld),
//        eMesh,
//        dataDir,
//        timeIndex,
//        nodeValues,
//        ensightCaseFile
//    );
//}


template<class Type>
Field<Type> map
(
    const Field<Type>& vf,
    const labelList& map1,
    const labelList& map2
)
{
    Field<Type> mf(map1.size() + map2.size());

    forAll(map1, i)
    {
        mf[i] = vf[map1[i]];
    }

    label offset = map1.size();

    forAll(map2, i)
    {
        mf[i + offset] = vf[map2[i]];
    }

    return mf;
}


template<class Type>
void writeField
(
    const char* key,
    const Field<Type>& vf,
    ensightStream& os
)
{
    if (returnReduce(vf.size(), sumOp<label>()) > 0)
    {
        if (Pstream::master())
        {
            os.write(key);

            for (direction i=0; i < pTraits<Type>::nComponents; ++i)
            {
                label cmpt = ensightPTraits<Type>::componentOrder[i];

                os.write(vf.component(cmpt));

                for (int slave=1; slave<Pstream::nProcs(); slave++)
                {
                    IPstream fromSlave(Pstream::scheduled, slave);
                    scalarField slaveData(fromSlave);
                    os.write(slaveData);
                }
            }
        }
        else
        {
            for (direction i=0; i < pTraits<Type>::nComponents; ++i)
            {
                label cmpt = ensightPTraits<Type>::componentOrder[i];

                OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
                toMaster<< vf.component(cmpt);
            }
        }
    }
}


template<class Type>
bool writePatchField
(
    const Field<Type>& pf,
    const label patchi,
    const label ensightPatchI,
    const faceSets& boundaryFaceSet,
    const ensightMesh::nFacePrimitives& nfp,
    ensightStream& os
)
{
    if (nfp.nTris || nfp.nQuads || nfp.nPolys)
    {
        if (Pstream::master())
        {
            os.writePartHeader(ensightPatchI);
        }

        writeField
        (
            "tria3",
            Field<Type>(pf, boundaryFaceSet.tris),
            os
        );

        writeField
        (
            "quad4",
            Field<Type>(pf, boundaryFaceSet.quads),
            os
        );

        writeField
        (
            "nsided",
            Field<Type>(pf, boundaryFaceSet.polys),
            os
        );

        return true;
    }
    else
    {
        return false;
    }
}


template<class Type>
void writePatchField
(
    const word& fieldName,
    const Field<Type>& pf,
    const word& patchName,
    const ensightMesh& eMesh,
    const fileName& dataDir,
    const label timeIndex,
    Ostream& ensightCaseFile
)
{
    const List<faceSets>& boundaryFaceSets = eMesh.boundaryFaceSets();
    const wordList& allPatchNames = eMesh.allPatchNames();
    const HashTable<ensightMesh::nFacePrimitives>&
        nPatchPrims = eMesh.nPatchPrims();

    label ensightPatchI = eMesh.patchPartOffset();

    label patchi = -1;

    forAll(allPatchNames, i)
    {
        if (allPatchNames[i] == patchName)
        {
            patchi = i;
            break;
        }
        ensightPatchI++;
    }

    ensightStream* filePtr(nullptr);
    if (Pstream::master())
    {
        // TODO: verify that these are indeed valid ensight variable names
        const word varName = patchName + '.' + fieldName;
        const fileName postFileName = ensightFile::subDir(timeIndex)/varName;

        // the data/ITER subdirectory must exist
        mkDir(dataDir/postFileName.path());

        if (timeIndex == 0)
        {
            const fileName dirName = dataDir.name()/ensightFile::mask();

            ensightCaseFile.setf(ios_base::left);
            ensightCaseFile
                << ensightPTraits<Type>::typeName << " per "
                << setw(20)
                << "element:"
                << " 1  "
                << setw(15)
                << varName.c_str() << ' '
                << (dirName/varName).c_str()
                << nl;
        }

        if (eMesh.format() == IOstream::BINARY)
        {
            filePtr = new ensightBinaryStream
            (
                dataDir/postFileName
            );
        }
        else
        {
            filePtr = new ensightAsciiStream
            (
                dataDir/postFileName
            );
        }

        filePtr->write(ensightPTraits<Type>::typeName);
    }

    ensightStream& os = *filePtr;

    if (patchi >= 0)
    {
        writePatchField
        (
            pf,
            patchi,
            ensightPatchI,
            boundaryFaceSets[patchi],
            nPatchPrims.find(patchName)(),
            os
        );
    }
    else
    {
        faceSets nullFaceSets;

        writePatchField
        (
            Field<Type>(),
            -1,
            ensightPatchI,
            nullFaceSets,
            nPatchPrims.find(patchName)(),
            os
        );
    }

    if (filePtr) // on master only
    {
        delete filePtr;
    }
}


template<class Type>
void ensightField
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const ensightMesh& eMesh,
    const fileName& dataDir,
    const label timeIndex,
    Ostream& ensightCaseFile
)
{
    Info<< ' ' << vf.name();

    const fvMesh& mesh = eMesh.mesh();

    const cellSets& meshCellSets = eMesh.meshCellSets();
    const List<faceSets>& boundaryFaceSets = eMesh.boundaryFaceSets();
    const wordList& allPatchNames = eMesh.allPatchNames();
    const wordHashSet& patchNames = eMesh.patchNames();
    const HashTable<ensightMesh::nFacePrimitives>&
        nPatchPrims = eMesh.nPatchPrims();
    const List<faceSets>& faceZoneFaceSets = eMesh.faceZoneFaceSets();
    const wordHashSet& faceZoneNames = eMesh.faceZoneNames();
    const HashTable<ensightMesh::nFacePrimitives>&
        nFaceZonePrims = eMesh.nFaceZonePrims();

    const labelList& tets = meshCellSets.tets;
    const labelList& pyrs = meshCellSets.pyrs;
    const labelList& prisms = meshCellSets.prisms;
    const labelList& wedges = meshCellSets.wedges;
    const labelList& hexes = meshCellSets.hexes;
    const labelList& polys = meshCellSets.polys;

    ensightStream* filePtr(nullptr);
    if (Pstream::master())
    {
        const ensight::VarName varName(vf.name());
        const fileName postFileName = ensightFile::subDir(timeIndex)/varName;

        // the data/ITER subdirectory must exist
        mkDir(dataDir/postFileName.path());

        if (timeIndex == 0)
        {
            const fileName dirName = dataDir.name()/ensightFile::mask();

            ensightCaseFile.setf(ios_base::left);
            ensightCaseFile
                << ensightPTraits<Type>::typeName
                << " per element:     1  "
                << setw(15)
                << varName.c_str() << ' '
                << (dirName/varName).c_str()
                << nl;
        }

        if (eMesh.format() == IOstream::BINARY)
        {
            filePtr = new ensightBinaryStream
            (
                dataDir/postFileName
            );
        }
        else
        {
            filePtr = new ensightAsciiStream
            (
                dataDir/postFileName
            );
        }

        filePtr->write(ensightPTraits<Type>::typeName);
    }

    ensightStream& os = *filePtr;

    if (patchNames.empty())
    {
        eMesh.barrier();

        if (Pstream::master())
        {
            os.writePartHeader(1);
        }

        writeField
        (
            "hexa8",
            map(vf, hexes, wedges),
            os
        );

        writeField
        (
            "penta6",
            Field<Type>(vf, prisms),
            os
        );

        writeField
        (
            "pyramid5",
            Field<Type>(vf, pyrs),
            os
        );

        writeField
        (
            "tetra4",
            Field<Type>(vf, tets),
            os
        );

        writeField
        (
            "nfaced",
            Field<Type>(vf, polys),
            os
        );
    }

    label ensightPatchI = eMesh.patchPartOffset();

    forAll(allPatchNames, patchi)
    {
        const word& patchName = allPatchNames[patchi];

        eMesh.barrier();

        if (patchNames.empty() || patchNames.found(patchName))
        {
            if
            (
                writePatchField
                (
                    vf.boundaryField()[patchi],
                    patchi,
                    ensightPatchI,
                    boundaryFaceSets[patchi],
                    nPatchPrims.find(patchName)(),
                    os
                )
            )
            {
                ensightPatchI++;
            }
        }
    }

    // write faceZones, if requested
    if (faceZoneNames.size())
    {
        // Interpolates cell values to faces - needed only when exporting
        // faceZones...
        GeometricField<Type, fvsPatchField, surfaceMesh> sf
        (
            linearInterpolate(vf)
        );

        forAllConstIter(wordHashSet, faceZoneNames, iter)
        {
            const word& faceZoneName = iter.key();

            eMesh.barrier();

            const label zoneID = mesh.faceZones().findZoneID(faceZoneName);
            const faceZone& fz = mesh.faceZones()[zoneID];

            // Prepare data to write
            label nIncluded = 0;
            forAll(fz, i)
            {
                if (eMesh.faceToBeIncluded(fz[i]))
                {
                    ++nIncluded;
                }
            }

            Field<Type> values(nIncluded);

            // Loop on the faceZone and store the needed field values
            label j = 0;
            forAll(fz, i)
            {
                label faceI = fz[i];
                if (mesh.isInternalFace(faceI))
                {
                    values[j] = sf[faceI];
                    ++j;
                }
                else
                {
                    if (eMesh.faceToBeIncluded(faceI))
                    {
                        label patchI = mesh.boundaryMesh().whichPatch(faceI);
                        const polyPatch& pp = mesh.boundaryMesh()[patchI];
                        label patchFaceI = pp.whichFace(faceI);
                        Type value = sf.boundaryField()[patchI][patchFaceI];
                        values[j] = value;
                        ++j;
                    }
                }
            }

            if
            (
                writePatchField
                (
                    values,
                    zoneID,
                    ensightPatchI,
                    faceZoneFaceSets[zoneID],
                    nFaceZonePrims.find(faceZoneName)(),
                    os
                )
            )
            {
                ensightPatchI++;
            }
        }
    }

    if (filePtr) // on master only
    {
        delete filePtr;
    }
}


template<class Type>
void ensightPointField
(
    const GeometricField<Type, pointPatchField, pointMesh>& pf,
    const ensightMesh& eMesh,
    const fileName& dataDir,
    const label timeIndex,
    Ostream& ensightCaseFile
)
{
    Info<< ' ' << pf.name();

    const fvMesh& mesh = eMesh.mesh();
    const wordList& allPatchNames = eMesh.allPatchNames();
    const wordHashSet& patchNames = eMesh.patchNames();
    const wordHashSet& faceZoneNames = eMesh.faceZoneNames();

    ensightStream* filePtr(nullptr);
    if (Pstream::master())
    {
        const ensight::VarName varName(pf.name());
        const fileName postFileName = ensightFile::subDir(timeIndex)/varName;

        // the data/ITER subdirectory must exist
        mkDir(dataDir/postFileName.path());

        if (timeIndex == 0)
        {
            const fileName dirName = dataDir.name()/ensightFile::mask();

            ensightCaseFile.setf(ios_base::left);
            ensightCaseFile
                << ensightPTraits<Type>::typeName
                << " per "
                << setw(20)
                << " node:"
                << " 1  "
                << setw(15)
                << varName.c_str() << ' '
                << (dirName/varName).c_str()
                << nl;
        }

        if (eMesh.format() == IOstream::BINARY)
        {
            filePtr = new ensightBinaryStream
            (
                dataDir/postFileName
            );
        }
        else
        {
            filePtr = new ensightAsciiStream
            (
                dataDir/postFileName
            );
        }

        filePtr->write(ensightPTraits<Type>::typeName);
    }

    ensightStream& os = *filePtr;

    if (eMesh.patchNames().empty())
    {
        eMesh.barrier();

        if (Pstream::master())
        {
            os.writePartHeader(1);
        }

        writeField
        (
            "coordinates",
            Field<Type>(pf.internalField(), eMesh.uniquePointMap()),
            os
        );
    }


    label ensightPatchI = eMesh.patchPartOffset();

    forAll(allPatchNames, patchi)
    {
        const word& patchName = allPatchNames[patchi];

        eMesh.barrier();

        if (patchNames.empty() || patchNames.found(patchName))
        {
            const fvPatch& p = mesh.boundary()[patchi];

            if (returnReduce(p.size(), sumOp<label>()) > 0)
            {
                // Renumber the patch points/faces into unique points
                labelList pointToGlobal;
                labelList uniqueMeshPointLabels;
                autoPtr<globalIndex> globalPointsPtr =
                mesh.globalData().mergePoints
                (
                    p.patch().meshPoints(),
                    p.patch().meshPointMap(),
                    pointToGlobal,
                    uniqueMeshPointLabels
                );

                if (Pstream::master())
                {
                    os.writePartHeader(ensightPatchI);
                }

                writeField
                (
                    "coordinates",
                    Field<Type>(pf.internalField(), uniqueMeshPointLabels),
                    os
                );

                ensightPatchI++;
            }
        }
    }

    // write faceZones, if requested
    if (faceZoneNames.size())
    {
        forAllConstIter(wordHashSet, faceZoneNames, iter)
        {
            const word& faceZoneName = iter.key();

            eMesh.barrier();

            const label zoneID = mesh.faceZones().findZoneID(faceZoneName);
            const faceZone& fz = mesh.faceZones()[zoneID];

            if (returnReduce(fz().nPoints(), sumOp<label>()) > 0)
            {
                // Renumber the faceZone points/faces into unique points
                labelList pointToGlobal;
                labelList uniqueMeshPointLabels;
                autoPtr<globalIndex> globalPointsPtr =
                mesh.globalData().mergePoints
                (
                    fz().meshPoints(),
                    fz().meshPointMap(),
                    pointToGlobal,
                    uniqueMeshPointLabels
                );

                if (Pstream::master())
                {
                    os.writePartHeader(ensightPatchI);
                }

                writeField
                (
                    "coordinates",
                    Field<Type>
                    (
                        pf.internalField(),
                        uniqueMeshPointLabels
                    ),
                    os
                );

                ensightPatchI++;
            }
        }
    }

    if (filePtr) // on master only
    {
        delete filePtr;
    }
}


template<class Type>
void ensightField
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const ensightMesh& eMesh,
    const fileName& dataDir,
    const label timeIndex,
    const bool nodeValues,
    Ostream& ensightCaseFile
)
{
    if (nodeValues)
    {
        tmp<GeometricField<Type, pointPatchField, pointMesh>> pfld
        (
            volPointInterpolation::New(vf.mesh()).interpolate(vf)
        );
        pfld.ref().rename(vf.name());

        ensightPointField<Type>
        (
            pfld,
            eMesh,
            dataDir,
            timeIndex,
            ensightCaseFile
        );
    }
    else
    {
        ensightField<Type>
        (
            vf,
            eMesh,
            dataDir,
            timeIndex,
            ensightCaseFile
        );
    }
}


// ************************************************************************* //
