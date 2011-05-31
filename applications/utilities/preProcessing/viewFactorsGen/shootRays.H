// All rays expressed as start face (local) index end end face (global)
// Pre-size by assuming a certain percentage is visible.

// Maximum lenght for dynamicList
const label maxDynListLenght = 10000;
//label lenghtCount = 0;

for (label procI = 0; procI < Pstream::nProcs(); procI++)
{
    // Shoot rays from me to procI. Note that even if processor has
    // 0 faces we still need to call findLine to keep calls synced.

    DynamicField<point> start(coarseMesh.nFaces());
    DynamicField<point> end(start.size());
    DynamicList<label> startIndex(start.size());
    DynamicList<label> endIndex(start.size());


    const pointField& myFc = remoteCoarseCf[Pstream::myProcNo()];
    const vectorField& myArea = remoteCoarseSf[Pstream::myProcNo()];

    const pointField& remoteArea = remoteCoarseSf[procI];
    const pointField& remoteFc = remoteCoarseCf[procI];

    if (myFc.size()*remoteFc.size() > 0)
    {
        forAll(myFc, i)
        {
            const point& fc = myFc[i];
            const vector& fA = myArea[i];

            forAll(remoteFc, j)
            {
                if (procI != Pstream::myProcNo() || i != j)
                {
                    const point& remFc = remoteFc[j];
                    const vector& remA = remoteArea[j];
                    const vector& d = remFc-fc;

                    if (((d & fA) < 0.) && ((d & remA) > 0))
                    {
                        //lenghtCount ++;
                        start.append(fc + 0.0001*d);
                        startIndex.append(i);
                        end.append(fc + 0.9999*d);
                        label globalI = globalNumbering.toGlobal(procI, j);
                        endIndex.append(globalI);

                        if (startIndex.size() > maxDynListLenght)
                        {
                            List<pointIndexHit> hitInfo(startIndex.size());
                            surfacesMesh.findLine
                            (
                                start,
                                end,
                                hitInfo
                            );
                            surfacesMesh.findLine(start, end, hitInfo);
                            forAll (hitInfo, rayI)
                            {
                                if (!hitInfo[rayI].hit())
                                {
                                    rayStartFace.append(startIndex[rayI]);
                                    rayEndFace.append(endIndex[rayI]);
                                }
                            }
                            //lenghtCount = 0;
                            start.clear();
                            startIndex.clear();
                            end.clear();
                            endIndex.clear();
                        }
                    }
                }
            }
        }
    }
    if (!start.empty())
    {
        List<pointIndexHit> hitInfo(startIndex.size());
        surfacesMesh.findLine
        (
            start,
            end,
            hitInfo
        );
        surfacesMesh.findLine(start, end, hitInfo);
        forAll (hitInfo, rayI)
        {
            if (!hitInfo[rayI].hit())
            {
                rayStartFace.append(startIndex[rayI]);
                rayEndFace.append(endIndex[rayI]);
            }
        }
    }
}
