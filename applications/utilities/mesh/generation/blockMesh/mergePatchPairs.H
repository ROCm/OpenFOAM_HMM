// Handle merging of patch pairs
{
    wordPairList mergePatchPairs;

    // Read in a list of dictionaries for the merge patch pairs
    if
    (
        meshDict.readIfPresent("mergePatchPairs", mergePatchPairs)
     && mergePatchPairs.size()
    )
    {
        Info<< "Creating merge patch pairs" << nl << endl;

        Info<< "Adding point and face zones" << endl;
        {
            auto& pzs = mesh.pointZones();
            pzs.clearAddressing();
            auto& fzs = mesh.faceZones();
            fzs.clearAddressing();

            forAll(mergePatchPairs, pairi)
            {
                const word mergeName
                (
                    mergePatchPairs[pairi].first()
                  + mergePatchPairs[pairi].second()
                  + name(pairi)
                );

                // An empty zone for cut points
                pzs.append
                (
                    new pointZone
                    (
                        mergeName + "CutPointZone",
                        pzs.size(),
                        pzs
                    )
                );

                // Master patch
                const word masterPatchName(mergePatchPairs[pairi].first());
                const polyPatch& masterPatch =
                    mesh.boundaryMesh()[masterPatchName];

                fzs.append
                (
                    new faceZone
                    (
                        mergeName + "MasterZone",
                        identity(masterPatch.range()),
                        false, // none are flipped
                        fzs.size(),
                        fzs
                    )
                );

                // Slave patch
                const word slavePatchName(mergePatchPairs[pairi].second());
                const polyPatch& slavePatch =
                    mesh.boundaryMesh()[slavePatchName];

                fzs.append
                (
                    new faceZone
                    (
                        mergeName + "SlaveZone",
                        identity(slavePatch.range()),
                        false, // none are flipped
                        fzs.size(),
                        fzs
                    )
                );

                // An empty zone for cut faces
                fzs.append
                (
                    new faceZone
                    (
                        mergeName + "CutFaceZone",
                        fzs.size(),
                        fzs
                    )
                );
            }  // end of all merge pairs
        }



        Info<< "Creating attachPolyTopoChanger" << endl;
        attachPolyTopoChanger polyMeshAttacher(mesh);
        polyMeshAttacher.setSize(mergePatchPairs.size());

        forAll(mergePatchPairs, pairi)
        {
            const word mergeName
            (
                mergePatchPairs[pairi].first()
              + mergePatchPairs[pairi].second()
              + name(pairi)
            );

            // Add the sliding interface mesh modifier
            polyMeshAttacher.set
            (
                pairi,
                new slidingInterface
                (
                    "couple" + name(pairi),
                    pairi,
                    polyMeshAttacher,
                    mergeName + "MasterZone",
                    mergeName + "SlaveZone",
                    mergeName + "CutPointZone",
                    mergeName + "CutFaceZone",
                    mergePatchPairs[pairi].first(),
                    mergePatchPairs[pairi].second(),
                    slidingInterface::INTEGRAL, // always integral
                    false,
                    intersection::VISIBLE
                )
            );
        }

        polyMeshAttacher.attach(true);
    }
    else
    {
        Info<< nl << "There are no merge patch pairs" << endl;
    }
}
