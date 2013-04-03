/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "GAMGAgglomeration.H"
#include "mapDistribute.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//template<class Type>
//void Foam::GAMGAgglomeration::gatherFaceField
//(
//    const label comm,
//    const labelList& procIDs,
//    const labelListList& procFaceMap,
//    const UList<Type>& fld,
//    List<Type>& allFld,
//    const int tag
//)
//{
//    if (Pstream::myProcNo(comm) == procIDs[0])
//    {
//        // Assign my local data
//        UIndirectList<Type>(allFld, procFaceMap[0]) = fld;
//
//        // Temporary storage
//        Field<Type> buff(allFld.size());
//
//        for (label i = 1; i < procIDs.size(); i++)
//        {
//            const labelList& map = procFaceMap[i];
//
//            IPstream::read
//            (
//                Pstream::scheduled,
//                procIDs[i],
//                reinterpret_cast<char*>(buff.begin()),
//                map.size()*sizeof(Type),
//                tag,
//                comm
//            );
//
//            forAll(map, i)
//            {
//                allFld[map[i]] = buff[i];
//            }
//        }
//    }
//    else
//    {
//        OPstream::write
//        (
//            Pstream::scheduled,
//            procIDs[0],
//            reinterpret_cast<const char*>(fld.begin()),
//            fld.byteSize(),
//            tag,
//            comm
//        );
//    }
//
//    //Pout<< "** GAMGAgglomeration: gatherFaceField :"
//    //    << " fld was:" << fld.size()
//    //    << " fld now:" << allFld.size() << endl;
//}
//
//
//template<class Type>
//void Foam::GAMGAgglomeration::gatherPatchFaceField
//(
//    const label comm,
//    const labelList& procIDs,
//    const labelListList& procBoundaryMap,
//    const labelListListList& procPatchFaceMap,
//
//    const List<List<Type> >& fld,
//    List<List<Type> >& allFld,
//    const int tag
//)
//{
//    if (Pstream::myProcNo(comm) == procIDs[0])
//    {
//        // Assign my local data
//        {
//            const labelList& patches = procBoundaryMap[0];
//
//            forAll(patches, intI)
//            {
//                label newIntI = patches[intI];
//                if (newIntI != -1)
//                {
//                    const labelList& map = procPatchFaceMap[0][intI];
//                    const List<Type>& patchFld = fld[intI];
//                    // Only copy unflipped faces. One of the processors will
//                    // have positive ones.
//                    forAll(map, i)
//                    {
//                        if (map[i] >= 0)
//                        {
//                            allFld[newIntI][map[i]] = patchFld[i];
//                        }
//                    }
//                }
//            }
//        }
//
//
//        // Temporary storage
//        List<List<Type> > buff(allFld.size());
//
//        for (label i = 1; i < procIDs.size(); i++)
//        {
//            IPstream fromSlave
//            (
//                Pstream::scheduled,
//                procIDs[i],
//                0,
//                tag,
//                comm
//            );
//            fromSlave >> buff;
//
//            const labelList& patches = procBoundaryMap[i];
//            forAll(patches, intI)
//            {
//                label newIntI = patches[intI];
//                if (newIntI != -1)
//                {
//                    const labelList& map = procPatchFaceMap[i][intI];
//                    const List<Type>& patchFld = buff[intI];
//                    forAll(map, i)
//                    {
//                        if (map[i] >= 0)
//                        {
//                            allFld[newIntI][map[i]] = patchFld[i];
//                        }
//                    }
//                }
//            }
//        }
//    }
//    else
//    {
//        OPstream toMaster
//        (
//            Pstream::scheduled,
//            procIDs[0],
//            0,
//            tag,
//            comm
//        );
//        toMaster << fld;
//    }
//
//    //Pout<< "** GAMGAgglomeration: gatherPatchFaceField :"
//    //    << " fld was:" << fld.size()
//    //    << " fld now:" << allFld.size() << endl;
//}


//XXXXXXXX

template<class Type>
void Foam::GAMGAgglomeration::mapField
(
    const faceMappingType mapType,
    const labelList& map,
    const UList<Type>& fld,
    UList<Type>& mappedFld
)
{
    switch (mapType)
    {
        case IGNORESIGN:
            forAll(map, i)
            {
                label allI = map[i];
                if (allI < 0)
                {
                    allI = -allI-1;
                }
                mappedFld[allI] = fld[i];
            }
        break;

        case APPLYSIGN:
            forAll(map, i)
            {
                label allI = map[i];
                if (allI < 0)
                {
                    allI = -allI-1;
                    mappedFld[allI] = -fld[i];
                }
                else
                {
                    mappedFld[allI] = fld[i];
                }
            }
        break;

        case DONTMAP:
            forAll(map, i)
            {
                label allI = map[i];
                if (allI >= 0)
                {
                    mappedFld[allI] = fld[i];
                }
            }
        break;

        default:
            FatalErrorIn("GAMGAgglomeration::mapField(..)")
                << "Unhandled mapping method." << exit(FatalError);
        break;
    }
}


template<class Type>
void Foam::GAMGAgglomeration::gatherFaceField
(
    const label comm,
    const labelList& procIDs,

    const faceMappingType mapType,

    // Field on internal faces
    const labelListList& procFaceMap,
    const UList<Type>& faceFld,

    // Field on patches faces
    const labelListList& procPatchMap,
    const labelListListList& procPatchFaceMap,
    const UList<List<Type> >& patchFld,

    List<Type>& allFaceFld,
    List<List<Type> >& allPatchFld,
    const int tag
)
{
    if (Pstream::myProcNo(comm) == procIDs[0])
    {
        // Map from local faces to faces
        mapField
        (
            mapType,
            procFaceMap[0],
            faceFld,
            allFaceFld
        );

        // Map patch fields
        const labelList& patches = procPatchMap[0];
        forAll(patches, intI)
        {
            const labelList& map = procPatchFaceMap[0][intI];

            label newIntI = patches[intI];
            if (newIntI != -1)
            {
                // Map into patch field
                mapField
                (
                    mapType,
                    map,
                    patchFld,
                    allPatchFld
                );
            }
            else
            {
                // Map into internal faces field
                mapField
                (
                    mapType,
                    map,
                    patchFld,
                    allFaceFld
                );
            }
        }



        // Temporary storage
        List<Type> procFaceFld;
        List<List<Type> > procPatchFld;

        for (label i = 1; i < procIDs.size(); i++)
        {
            IPstream fromSlave
            (
                Pstream::scheduled,
                procIDs[i],
                0,
                tag,
                comm
            );
            fromSlave >> procFaceFld >> procPatchFld;

            // Map from faces to faces
            mapField
            (
                mapType,
                procFaceMap[i],
                procFaceFld,
                allFaceFld
            );


            const labelList& patches = procPatchMap[i];
            forAll(patches, intI)
            {
                const List<Type>& patchFld = procPatchFld[intI];
                const labelList& map = procPatchFaceMap[i][intI];

                label newIntI = patches[intI];
                if (newIntI != -1)
                {
                    // Map into patch field
                    mapField
                    (
                        mapType,
                        map,
                        patchFld,
                        allPatchFld[newIntI]
                    );
                }
                else
                {
                    // Map into internal faces field
                    mapField
                    (
                        mapType,
                        map,
                        patchFld,
                        allFaceFld
                    );
                }
            }
        }
    }
    else
    {
        OPstream toMaster
        (
            Pstream::scheduled,
            procIDs[0],
            0,
            tag,
            comm
        );
        toMaster << faceFld << patchFld;
    }
}
template<class Type>
void Foam::GAMGAgglomeration::gatherList
(
    const label comm,
    const labelList& procIDs,

    const Type& myVal,
    List<Type>& allVals,
    const int tag
)
{
    if (Pstream::myProcNo(comm) == procIDs[0])
    {
        allVals.setSize(procIDs.size());

        allVals[0] = myVal;
        for (label i = 1; i < procIDs.size(); i++)
        {
            IPstream fromSlave
            (
                Pstream::scheduled,
                procIDs[i],
                0,
                tag,
                comm
            );

            fromSlave >> allVals[i];
        }
    }
    else
    {
        OPstream toMaster
        (
            Pstream::scheduled,
            procIDs[0],
            0,
            tag,
            comm
        );
        toMaster << myVal;
    }
}


template<class Type>
void Foam::GAMGAgglomeration::restrictField
(
    Field<Type>& cf,
    const Field<Type>& ff,
    const labelList& fineToCoarse
) const
{
    cf = pTraits<Type>::zero;

    forAll(ff, i)
    {
        cf[fineToCoarse[i]] += ff[i];
    }
}
//XXXXXXXX




template<class Type>
void Foam::GAMGAgglomeration::restrictField
(
    Field<Type>& cf,
    const Field<Type>& ff,
    const label fineLevelIndex,
    const bool procAgglom
) const
{
    const labelList& fineToCoarse = restrictAddressing_[fineLevelIndex];

    if (!procAgglom && ff.size() != fineToCoarse.size())
    {
        FatalErrorIn
        (
            "void GAMGAgglomeration::restrictField"
            "(Field<Type>& cf, const Field<Type>& ff, "
            "const label fineLevelIndex) const"
        )   << "field does not correspond to level " << fineLevelIndex
            << " sizes: field = " << ff.size()
            << " level = " << fineToCoarse.size()
            << abort(FatalError);
    }

    Pout<< "** GAMGAgglomeration: restrictField :"
        << " fineLevel:" << fineLevelIndex
        << " agglomerating from:" << ff.size()
        << " down to:" << cf.size()
        << " max fineToCoarse:" << max(fineToCoarse) << endl;

    cf = pTraits<Type>::zero;

    forAll(ff, i)
    {
        cf[fineToCoarse[i]] += ff[i];
    }


    label coarseLevelIndex = fineLevelIndex+1;

    if (procAgglom && hasProcMesh(coarseLevelIndex))
    {
        label fineComm = UPstream::parent(procCommunicator_[coarseLevelIndex]);

        Pout<< "** GAMGAgglomeration: restrictField :"
            << " fineLevel:" << fineLevelIndex
            << " coarseLevel:" << coarseLevelIndex
            << " have proc agglomeration with comm:"
            << fineComm << endl;


        const List<int>& procIDs = agglomProcIDs(coarseLevelIndex);
        const labelList& offsets = cellOffsets(coarseLevelIndex);

        Pout<< "procIDs:" << procIDs << endl;
        Pout<< "offsets:" << offsets << endl;
        Pout<< "myProcNo:" << UPstream::myProcNo(fineComm) << endl;

        const globalIndex cellOffsetter(offsets);

        Pout<< "** GAMGAgglomeration: restrictField :"
            << " cf was:" << cf.size();

        cellOffsetter.gather(fineComm, procIDs, cf);

        Pout<< " cf now:" << cf.size() << endl;
    }
}


template<class Type>
void Foam::GAMGAgglomeration::restrictFaceField
(
    Field<Type>& cf,
    const Field<Type>& ff,
    const label fineLevelIndex
) const
{
    const labelList& fineToCoarse = faceRestrictAddressing_[fineLevelIndex];

    if (ff.size() != fineToCoarse.size())
    {
        FatalErrorIn
        (
            "void GAMGAgglomeration::restrictFaceField"
            "(Field<Type>& cf, const Field<Type>& ff, "
            "const label fineLevelIndex) const"
        )   << "field does not correspond to level " << fineLevelIndex
            << " sizes: field = " << ff.size()
            << " level = " << fineToCoarse.size()
            << abort(FatalError);
    }

    cf = pTraits<Type>::zero;

    forAll(fineToCoarse, ffacei)
    {
        label cFace = fineToCoarse[ffacei];

        if (cFace >= 0)
        {
            cf[cFace] += ff[ffacei];
        }
    }
}


//template<class Type>
//void Foam::GAMGAgglomeration::restrictFaceField
//(
//    Field<Type>& cf,
//    const Field<Type>& ff,
//    const label fineLevelIndex,
//    const bool dummy
//) const
//{
//    restrictFaceField(cf, ff, fineLevelIndex);
//
//    if (hasProcMesh(fineLevelIndex))
//    {
//        label fineComm = UPstream::parent(procCommunicator_[fineLevelIndex]);
//
//        Pout<< "** GAMGAgglomeration: restrictFaceField :"
//            << " fineLevel:" << fineLevelIndex
//            << " coarseLevel:" << fineLevelIndex+1
//            << " have proc agglomeration with comm:"
//            << fineComm << endl;
//
//        const labelListList& map = faceMap(fineLevelIndex);
//        const List<int>& procIDs = agglomProcIDs(fineLevelIndex);
//
//        Field<Type> allCf;
//        if (Pstream::myProcNo(fineComm) == procIDs[0])
//        {
//            allCf.setSize
//            (
//                meshLevel(fineLevelIndex).lduAddr().size()
//            );
//        }
//
//        gatherFaceField(fineComm, procIDs, map, cf, allCf);
//
//        Pout<< "** GAMGAgglomeration: restrictFaceField :"
//            << " cf was:" << cf.size()
//            << " cf now:" << allCf.size() << endl;
//
//        cf.transfer(allCf);
//    }
//}


template<class Type>
void Foam::GAMGAgglomeration::prolongField
(
    Field<Type>& ff,
    const Field<Type>& cf,
    const label levelIndex,
    const bool procAgglom
) const
{
    const labelList& fineToCoarse = restrictAddressing_[levelIndex];

    Pout<< "** GAMGAgglomeration: prolongField :"
        << " fineLevel:" << levelIndex
        << " prolonging from:" << cf.size()
        << " to:" << ff.size()
        << " max fineToCoarse:" << max(fineToCoarse) << endl;

    label coarseLevelIndex = levelIndex+1;

    if (procAgglom && hasProcMesh(coarseLevelIndex))
    {
        label coarseComm = UPstream::parent
        (
            procCommunicator_[coarseLevelIndex]
        );

        Pout<< "** GAMGAgglomeration: prolongField :"
            << " fineLevel:" << levelIndex
            << " coarseLevel:" << coarseLevelIndex
            << " have proc agglomeration with comm:"
            << coarseComm << endl;

        const List<int>& procIDs = agglomProcIDs(coarseLevelIndex);
        const labelList& offsets = cellOffsets(coarseLevelIndex);

        Pout<< "procIDs:" << procIDs << endl;
        Pout<< "offsets:" << offsets << endl;
        Pout<< "myProcNo:" << UPstream::myProcNo(coarseComm) << endl;

        const globalIndex cellOffsetter(offsets);

        label localSize = nCells_[levelIndex];
        Pout<< "localSize:" << localSize << endl;


        Field<Type> allCf(localSize);
        cellOffsetter.scatter
        (
            coarseComm,
            procIDs,
            cf,
            allCf
        );

        Pout<< "** GAMGAgglomeration: prolongField :"
            << " cf was:" << cf.size()
            << " cf now:" << allCf.size() << endl;

        forAll(fineToCoarse, i)
        {
            ff[i] = allCf[fineToCoarse[i]];
        }
    }
    else
    {
        forAll(fineToCoarse, i)
        {
            ff[i] = cf[fineToCoarse[i]];
        }
    }
}


//template<class Type>
//void Foam::GAMGAgglomeration::scatter
//(
//    Field<Type>& ff,
//    const Field<Type>& allFf,
//    const label levelIndex
//) const
//{
//    label coarseComm = UPstream::parent
//    (
//        procCommunicator_[levelIndex]
//    );
//
//    Pout<< "** GAMGAgglomeration: prolongField :"
//        << " fineLevel:" << levelIndex
//        << " have proc agglomeration with comm:"
//        << coarseComm << endl;
//
//    const List<int>& procIDs = agglomProcIDs(levelIndex);
//    const labelList& offsets = cellOffsets(levelIndex);
//
//    Pout<< "procIDs:" << procIDs << endl;
//    Pout<< "offsets:" << offsets << endl;
//    Pout<< "myProcNo:" << UPstream::myProcNo(coarseComm) << endl;
//
//    const globalIndex cellOffsetter(offsets);
//
//    label localSize = nCells_[levelIndex-1];
//    Pout<< "localSize:" << localSize << endl;
//
//
//    ff.setSize(localSize);
//    cellOffsetter.scatter
//    (
//        coarseComm,
//        procIDs,
//        allFf,
//        ff
//    );
//
//    Pout<< "** GAMGAgglomeration: prolongField :"
//        << " cf was:" << allFf.size()
//        << " cf now:" << ff.size() << endl;
//}



// ************************************************************************* //
