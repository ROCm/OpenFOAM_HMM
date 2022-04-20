// All rays expressed as start face (local) index end end face (global)
// Pre-size by assuming a certain percentage is visible.

// Maximum length for dynamicList
const label maxDynListLength
(
    viewFactorDict.getOrDefault<label>("maxDynListLength", 1000000000)
);

for (const int proci : Pstream::allProcs())
{

    std::vector<pbrt::Point3f> start;
    start.reserve(coarseMesh.nFaces());
    std::vector<pbrt::Vector3f> dir;
    dir.reserve(start.size());

    std::vector<pbrt::Point3f> end(start.size());
    end.reserve(start.size());

    DynamicList<label> startIndex(start.size());
    DynamicList<label> endIndex(start.size());

    DynamicList<label> startAgg(start.size());
    DynamicList<label> endAgg(start.size());

    const pointField& myFc = remoteCoarseCf[Pstream::myProcNo()];
    const vectorField& myArea = remoteCoarseSf[Pstream::myProcNo()];
    const labelField& myAgg = remoteCoarseAgg[Pstream::myProcNo()];

    const pointField& remoteArea = remoteCoarseSf[proci];
    const pointField& remoteFc = remoteCoarseCf[proci];
    const labelField& remoteAgg = remoteCoarseAgg[proci];

    pbrt::Vector3f smallV(SMALL, SMALL, SMALL);

    label i = 0;
    label j = 0;
    do
    {
        for (; i < myFc.size(); i++)
        {
            const point& fc = myFc[i];
            const vector& fA = myArea[i];
            const label& fAgg = myAgg[i];

            for (; j < remoteFc.size(); j++)
            {
                if (proci != Pstream::myProcNo() || i != j)
                {
                    const point& remFc = remoteFc[j];
                    const vector& remA = remoteArea[j];
                    const label& remAgg = remoteAgg[j];

                    const vector d = remFc - fc;

                    const vector nd = d/mag(d);
                    const vector nfA = fA/mag(fA);
                    const vector nremA = remA/mag(remA);

                    if (((nd & nfA) < 0) && ((nd & nremA) > 0))
                    {
                        pbrt::Vector3f direction(d[0], d[1], d[2]);
                        direction += smallV;
                        pbrt::Point3f s(fc[0], fc[1], fc[2]);
                        end.push_back(s + direction);

                        s += pbrt::Point3f(0.01*d[0], 0.01*d[1], 0.01*d[2]);

                        start.push_back(s);
                        startIndex.append(i);
                        startAgg.append(globalNumbering.toGlobal(proci, fAgg));

                        dir.push_back(direction);

                        label globalI = globalNumbering.toGlobal(proci, j);
                        endIndex.append(globalI);
                        endAgg.append(globalNumbering.toGlobal(proci, remAgg));
                        if (startIndex.size() > maxDynListLength)
                        {
                            FatalErrorInFunction
                                << "Dynamic list need from capacity."
                                << "Actual size maxDynListLength : "
                                <<  maxDynListLength
                                << abort(FatalError);
                        }
                    }
                }
            }

            if (j == remoteFc.size())
            {
                j = 0;
            }
        }

    } while (i < myFc.size());

    label totalRays(startIndex.size());
    reduce(totalRays, sumOp<scalar>());

    if (Pstream::master() && debug)
    {
        Info<< "Number of rays :" << totalRays << endl;
    }

    DynamicList<label> dRayIs;

    std::vector<pbrt::Ray> raysInt;
    raysInt.reserve(start.size());

    std::vector<int> raysTriIndex;
    raysTriIndex.reserve(start.size());

    for (unsigned long rayI = 0; rayI < start.size(); ++rayI)
    {
        pbrt::Ray ray(start[rayI], dir[rayI]);

        pbrt::SurfaceInteraction isect;
        bool intersects = accel->Intersect(ray, &isect);
        const Triangle* tri = dynamic_cast<const Triangle *>(isect.shape);

        if (intersects)
        {
            const int index = tri->faceIndex;

            const Vector3f delta(ray(ray.tMax) - end[rayI]);

            if (delta.Length() < intTol)
            {
                rayStartFace.append(startIndex[rayI]);
                rayEndFace.append(endIndex[rayI]);
            }
            else if (triSurfaceToAgglom[index] == startAgg[rayI])
            {
                dRayIs.append(rayI);
                raysInt.push_back(ray);
                raysTriIndex.push_back(index);
            }
        }
    }
    start.clear();

    labelField rayIs;
    rayIs.transfer(dRayIs);
    dRayIs.clear();

    boolList converged(rayIs.size(), false);
    label nConverged = 0;
    label iter = 0;
    do
    {
        forAll(rayIs, rayI)
        {
            const label rayID = rayIs[rayI];
            pbrt::Ray& rayHit = raysInt[rayI];
            const int index =  raysTriIndex[rayI];

            if (!converged[rayI])
            {
                if (triSurfaceToAgglom[index] == startAgg[rayID])
                {
                    rayHit.tMax *= 1.1;
                    rayHit.o = rayHit(rayHit.tMax);

                    pbrt::SurfaceInteraction isect;
                    bool intersects = accel->Intersect(rayHit, &isect);
                    if (intersects)
                    {
                        const Triangle* tri =
                            dynamic_cast<const Triangle *>(isect.shape);
                        const int index = tri->faceIndex;

                        raysTriIndex[rayI] = index;
                    }
                }
                else if (triSurfaceToAgglom[index] == endAgg[rayID])
                {
                    converged[rayI] = true;
                    nConverged++;
                    rayStartFace.append(startIndex[rayID]);
                    rayEndFace.append(endIndex[rayID]);
                }
            }
        }

        iter ++;

    } while (nConverged < converged.size() && iter < 10);

    startIndex.clear();
    end.clear();
    endIndex.clear();
    startAgg.clear();
    endAgg.clear();
}
