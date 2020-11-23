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

        // Create and add point and face zones and mesh modifiers
        List<pointZone*> pz(mergePatchPairs.size());
        List<faceZone*> fz(3*mergePatchPairs.size());
        List<cellZone*> cz;

        forAll(mergePatchPairs, pairi)
        {
            const word mergeName
            (
                mergePatchPairs[pairi].first()
              + mergePatchPairs[pairi].second()
              + name(pairi)
            );

            // An empty zone for cut points
            pz[pairi] = new pointZone
            (
                mergeName + "CutPointZone",
                0,
                mesh.pointZones()
            );

            // Master patch
            const word masterPatchName(mergePatchPairs[pairi].first());
            const polyPatch& masterPatch =
                mesh.boundaryMesh()[masterPatchName];

            fz[3*pairi] = new faceZone
            (
                mergeName + "MasterZone",
                identity(masterPatch.range()),
                false, // none are flipped
                0,
                mesh.faceZones()
            );

            // Slave patch
            const word slavePatchName(mergePatchPairs[pairi].second());
            const polyPatch& slavePatch =
                mesh.boundaryMesh()[slavePatchName];

            fz[3*pairi + 1] = new faceZone
            (
                mergeName + "SlaveZone",
                identity(slavePatch.range()),
                false, // none are flipped
                1,
                mesh.faceZones()
            );

            // An empty zone for cut faces
            fz[3*pairi + 2] = new faceZone
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
