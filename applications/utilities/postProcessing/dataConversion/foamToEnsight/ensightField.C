/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "ensightField.H"
#include "fvMesh.H"
#include "volFields.H"
#include "OFstream.H"
#include "IOmanip.H"
#include "itoa.H"
#include "volPointInterpolation.H"
#include "ensightBinaryStream.H"
#include "ensightAsciiStream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
scalarField map
(
    const Field<Type>& vf,
    const labelList& map,
    const label cmpt
)
{
    scalarField mf(map.size());

    forAll(map, i)
    {
        mf[i] = component(vf[map[i]], cmpt);
    }

    return mf;
}


template<class Type>
scalarField map
(
    const Field<Type>& vf,
    const labelList& map1,
    const labelList& map2,
    const label cmpt
)
{
    scalarField mf(map1.size() + map2.size());

    forAll(map1, i)
    {
        mf[i] = component(vf[map1[i]], cmpt);
    }

    label offset = map1.size();

    forAll(map2, i)
    {
        mf[i + offset] = component(vf[map2[i]], cmpt);
    }

    return mf;
}


template<class Type>
void writeAllData
(
    const char* key,
    const Field<Type>& vf,
    const labelList& prims,
    const label nPrims,
    ensightStream& ensightFile
)
{
    if (nPrims)
    {
        if (Pstream::master())
        {
            ensightFile.write(key);

            for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
            {
                scalarField masterData(map(vf, prims, cmpt));
                ensightFile.write(masterData);

                for (int slave=1; slave<Pstream::nProcs(); slave++)
                {
                    IPstream fromSlave(Pstream::scheduled, slave);
                    scalarField slaveData(fromSlave);
                    ensightFile.write(slaveData);
                }
            }
        }
        else
        {
            for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
            {
                OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
                toMaster<< map(vf, prims, cmpt);
            }
        }
    }
}


template<class Type>
void writeAllFaceData
(
    const char* key,
    const labelList& prims,
    const label nPrims,
    const Field<Type>& pf,
    ensightStream& ensightFile
)
{
    if (nPrims)
    {
        if (Pstream::master())
        {
            ensightFile.write(key);

            for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
            {
                ensightFile.write(map(pf, prims, cmpt));

                for (int slave=1; slave<Pstream::nProcs(); slave++)
                {
                    IPstream fromSlave(Pstream::scheduled, slave);
                    scalarField pf(fromSlave);

                    ensightFile.write(pf);
                }
            }
        }
        else
        {
            for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
            {
                OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
                toMaster<< map(pf, prims, cmpt);
            }
        }
    }
}


template<class Type>
bool writePatchField
(
    const Foam::Field<Type>& pf,
    const Foam::label patchi,
    const Foam::label ensightPatchI,
    const Foam::faceSets& boundaryFaceSet,
    const Foam::ensightMesh::nFacePrimitives& nfp,
    ensightStream& ensightFile
)
{
    if (nfp.nTris || nfp.nQuads || nfp.nPolys)
    {
        if (Pstream::master())
        {
            ensightFile.writePartHeader(ensightPatchI);
        }

        writeAllFaceData
        (
            "tria3",
            boundaryFaceSet.tris,
            nfp.nTris,
            pf,
            ensightFile
        );

        writeAllFaceData
        (
            "quad4",
            boundaryFaceSet.quads,
            nfp.nQuads,
            pf,
            ensightFile
        );

        writeAllFaceData
        (
            "nsided",
            boundaryFaceSet.polys,
            nfp.nPolys,
            pf,
            ensightFile
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
    const Foam::word& fieldName,
    const Foam::Field<Type>& pf,
    const Foam::word& patchName,
    const Foam::ensightMesh& eMesh,
    const Foam::fileName& postProcPath,
    const Foam::word& prepend,
    const Foam::label timeIndex,
    const bool binary,
    Foam::Ostream& ensightCaseFile
)
{
    const Time& runTime = eMesh.mesh().time();

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


    word pfName = patchName + '.' + fieldName;

    word timeFile = prepend + itoa(timeIndex);

    autoPtr<ensightStream> ensightFilePtr;
    if (Pstream::master())
    {
        if (timeIndex == 0)
        {
            ensightCaseFile.setf(ios_base::left);

            ensightCaseFile
                << pTraits<Type>::typeName
                << " per element:            1       "
                << setw(15) << pfName
                << (' ' + prepend + "***." + pfName).c_str()
                << nl;
        }

        // set the filename of the ensight file
        fileName ensightFileName(timeFile + "." + pfName);

        if (binary)
        {
            ensightFilePtr.reset
            (
                new ensightBinaryStream
                (
                    postProcPath/ensightFileName,
                    runTime
                )
            );
        }
        else
        {
            ensightFilePtr.reset
            (
                new ensightAsciiStream
                (
                    postProcPath/ensightFileName,
                    runTime
                )
            );
        }
    }

    ensightStream& ensightFile = ensightFilePtr();

    if (Pstream::master())
    {
        ensightFile.write(pTraits<Type>::typeName);
    }

    if (patchi >= 0)
    {
        writePatchField
        (
            pf,
            patchi,
            ensightPatchI,
            boundaryFaceSets[patchi],
            nPatchPrims.find(patchName)(),
            ensightFile
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
            ensightFile
        );
    }
}


template<class Type>
void ensightField
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const Foam::ensightMesh& eMesh,
    const Foam::fileName& postProcPath,
    const Foam::word& prepend,
    const Foam::label timeIndex,
    const bool binary,
    Foam::Ostream& ensightCaseFile
)
{
    Info<< "Converting field " << vf.name() << endl;

    word timeFile = prepend + itoa(timeIndex);

    const fvMesh& mesh = eMesh.mesh();
    const Time& runTime = mesh.time();

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

    autoPtr<ensightStream> ensightFilePtr;
    if (Pstream::master())
    {
        // set the filename of the ensight file
        fileName ensightFileName(timeFile + "." + vf.name());

        if (binary)
        {
            ensightFilePtr.reset
            (
                new ensightBinaryStream
                (
                    postProcPath/ensightFileName,
                    runTime
                )
            );
        }
        else
        {
            ensightFilePtr.reset
            (
                new ensightAsciiStream
                (
                    postProcPath/ensightFileName,
                    runTime
                )
            );
        }
    }

    ensightStream& ensightFile = ensightFilePtr();

    if (patchNames.empty())
    {
        eMesh.barrier();

        if (Pstream::master())
        {
            if (timeIndex == 0)
            {
                ensightCaseFile.setf(ios_base::left);

                ensightCaseFile
                    << pTraits<Type>::typeName
                    << " per element:            1       "
                    << setw(15) << vf.name()
                    << (' ' + prepend + "***." + vf.name()).c_str()
                    << nl;
            }

            ensightFile.write(pTraits<Type>::typeName);
            ensightFile.writePartHeader(1);
        }

        if (meshCellSets.nHexesWedges)
        {
            if (Pstream::master())
            {
                ensightFile.write("hexa8");

                for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
                {
                    scalarField masterData(map(vf, hexes, wedges, cmpt));
                    ensightFile.write(masterData);

                    for (int slave=1; slave<Pstream::nProcs(); slave++)
                    {
                        IPstream fromSlave(Pstream::scheduled, slave);
                        scalarField data(fromSlave);
                        ensightFile.write(data);
                    }
                }
            }
            else
            {
                for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
                {
                    OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
                    toMaster<< map(vf, hexes, wedges, cmpt);
                }
            }
        }

        writeAllData
        (
            "penta6",
            vf,
            prisms,
            meshCellSets.nPrisms,
            ensightFile
        );

        writeAllData
        (
            "pyramid5",
            vf,
            pyrs,
            meshCellSets.nPyrs,
            ensightFile
        );

        writeAllData
        (
            "tetra4",
            vf,
            tets,
            meshCellSets.nTets,
            ensightFile
        );

        writeAllData
        (
            "nfaced",
            vf,
            polys,
            meshCellSets.nPolys,
            ensightFile
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
                    ensightFile
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

            label zoneID = mesh.faceZones().findZoneID(faceZoneName);

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
                    ensightFile
                )
            )
            {
                ensightPatchI++;
            }
        }
    }
}


template<class Type>
void ensightPointField
(
    const GeometricField<Type, pointPatchField, pointMesh>& pf,
    const Foam::ensightMesh& eMesh,
    const Foam::fileName& postProcPath,
    const Foam::word& prepend,
    const Foam::label timeIndex,
    const bool binary,
    Foam::Ostream& ensightCaseFile
)
{
    Info<< "Converting field " << pf.name() << endl;

    word timeFile = prepend + itoa(timeIndex);

    //const fvMesh& mesh = eMesh.mesh();
    //const Time& runTime = mesh.time();

    autoPtr<ensightStream> ensightFilePtr;
    if (Pstream::master())
    {
        // set the filename of the ensight file
        fileName ensightFileName(timeFile + "." + pf.name());

        if (binary)
        {
            ensightFilePtr.reset
            (
                new ensightBinaryStream
                (
                    postProcPath/ensightFileName,
                    eMesh.mesh().time()
                )
            );
        }
        else
        {
            ensightFilePtr.reset
            (
                new ensightAsciiStream
                (
                    postProcPath/ensightFileName,
                    eMesh.mesh().time()
                )
            );
        }
    }

    ensightStream& ensightFile = ensightFilePtr();

    if (eMesh.patchNames().empty())
    {
        eMesh.barrier();

        if (Pstream::master())
        {
            if (timeIndex == 0)
            {
                ensightCaseFile.setf(ios_base::left);

                ensightCaseFile
                    << pTraits<Type>::typeName
                    << " per node:            1       "
                    << setw(15) << pf.name()
                    << (' ' + prepend + "***." + pf.name()).c_str()
                    << nl;
            }

            ensightFile.write(pTraits<Type>::typeName);
            ensightFile.write("part");
            ensightFile.write(1);
        }

        if (Pstream::master())
        {
            ensightFile.write("coordinates");

            Field<Type> uniqueFld(pf.internalField(), eMesh.uniquePointMap());

            for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
            {
                ensightFile.write(uniqueFld.component(cmpt));

                for (int slave=1; slave<Pstream::nProcs(); slave++)
                {
                    IPstream fromSlave(Pstream::scheduled, slave);
                    scalarField data(fromSlave);
                    ensightFile.write(data);
                }
            }
        }
        else
        {
            Field<Type> uniqueFld(pf.internalField(), eMesh.uniquePointMap());
            for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
            {
                OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
                toMaster<< uniqueFld.component(cmpt);
            }
        }
    }
}


template<class Type>
void ensightField
(
    const Foam::IOobject& fieldObject,
    const Foam::ensightMesh& eMesh,
    const Foam::fileName& postProcPath,
    const Foam::word& prepend,
    const Foam::label timeIndex,
    const bool binary,
    const bool nodeValues,
    Foam::Ostream& ensightCaseFile
)
{
    // Read field
    GeometricField<Type, fvPatchField, volMesh> vf(fieldObject, eMesh.mesh());

    if (nodeValues)
    {
        tmp<GeometricField<Type, pointPatchField, pointMesh> > pfld
        (
            volPointInterpolation::New(eMesh.mesh()).interpolate(vf)
        );
        pfld().rename(vf.name());

        ensightPointField<Type>
        (
            pfld,
            eMesh,
            postProcPath,
            prepend,
            timeIndex,
            binary,
            ensightCaseFile
        );
    }
    else
    {
        ensightField<Type>
        (
            vf,
            eMesh,
            postProcPath,
            prepend,
            timeIndex,
            binary,
            ensightCaseFile
        );
    }
}


// ************************************************************************* //
