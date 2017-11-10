        if (mergePatchPairs.size())
        {
            Info<< "Creating merge patch pairs" << nl << endl;

            // Create and add point and face zones and mesh modifiers
            List<pointZone*> pz(mergePatchPairs.size());
            List<faceZone*> fz(3*mergePatchPairs.size());
            List<cellZone*> cz(0);

            forAll(mergePatchPairs, pairI)
            {
                const word mergeName
                (
                    mergePatchPairs[pairI].first()
                  + mergePatchPairs[pairI].second()
                  + name(pairI)
                );

                // An empty zone for cut points
                pz[pairI] = new pointZone
                (
                    mergeName + "CutPointZone",
                    0,
                    mesh.pointZones()
                );

                // Master patch
                const word masterPatchName(mergePatchPairs[pairI].first());
                const polyPatch& masterPatch =
                    mesh.boundaryMesh()[masterPatchName];

                labelList isf(masterPatch.size());

                forAll(isf, i)
                {
                    isf[i] = masterPatch.start() + i;
                }

                fz[3*pairI] = new faceZone
                (
                    mergeName + "MasterZone",
                    isf,
                    false, // none are flipped
                    0,
                    mesh.faceZones()
                );

                // Slave patch
                const word slavePatchName(mergePatchPairs[pairI].second());
                const polyPatch& slavePatch =
                    mesh.boundaryMesh()[slavePatchName];

                labelList osf(slavePatch.size());

                forAll(osf, i)
                {
                    osf[i] = slavePatch.start() + i;
                }

                fz[3*pairI + 1] = new faceZone
                (
                    mergeName + "SlaveZone",
                    osf,
                    false, // none are flipped
                    1,
                    mesh.faceZones()
                );

                // An empty zone for cut faces
                fz[3*pairI + 2] = new faceZone
                (
                    mergeName + "CutFaceZone",
                    2,
                    mesh.faceZones()
                );
            }  // end of all merge pairs

            Info<< "Adding point and face zones" << endl;
            mesh.addZones(pz, fz, cz);


            Info<< "Creating attachPolyTopoChanger" << endl;
            attachPolyTopoChanger polyMeshAttacher(mesh);
            polyMeshAttacher.setSize(mergePatchPairs.size());

            forAll(mergePatchPairs, pairI)
            {
                const word mergeName
                (
                    mergePatchPairs[pairI].first()
                  + mergePatchPairs[pairI].second()
                  + name(pairI)
                );

                // Add the sliding interface mesh modifier
                polyMeshAttacher.set
                (
                    pairI,
                    new slidingInterface
                    (
                        "couple" + name(pairI),
                        pairI,
                        polyMeshAttacher,
                        mergeName + "MasterZone",
                        mergeName + "SlaveZone",
                        mergeName + "CutPointZone",
                        mergeName + "CutFaceZone",
                        mergePatchPairs[pairI].first(),
                        mergePatchPairs[pairI].second(),
                        slidingInterface::INTEGRAL, // always integral
                        intersection::VISIBLE
                    )
                );
            }

            polyMeshAttacher.attach(true);
        }
