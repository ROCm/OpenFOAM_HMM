// Set any cellZones
// Note cell labelling unaffected by previous mergePatchPairs

{
    const label nZones = blocks.numZonedBlocks();
    if (nZones)
    {
        Info<< nl << "Adding cell zones" << endl;

        // Map from zoneName to cellZone index
        HashTable<label> zoneMap(nZones);

        // Cells per zone.
        List<DynamicList<label>> zoneCells(nZones);

        // Running cell counter
        label celli = 0;

        // Largest zone so far
        label freeZoneI = 0;

        for (const block& b : blocks)
        {
            const word& zoneName = b.zoneName();
            const label nCellsInBlock = b.cells().size();

            if (zoneName.size())
            {
                const auto iter = zoneMap.cfind(zoneName);

                label zonei;

                if (iter.found())
                {
                    zonei = *iter;
                }
                else
                {
                    zonei = freeZoneI++;

                    Info<< "    " << zonei << '\t' << zoneName << endl;

                    zoneMap.insert(zoneName, zonei);
                }


                // Fill with cell ids

                zoneCells[zonei].reserve
                (
                    zoneCells[zonei].size() + nCellsInBlock
                );

                const label endOfFill = celli + nCellsInBlock;

                for (; celli < endOfFill; ++celli)
                {
                    zoneCells[zonei].append(celli);
                }
            }
            else
            {
                celli += nCellsInBlock;
            }
        }

        List<cellZone*> cz(zoneMap.size());
        forAllConstIters(zoneMap, iter)
        {
            const word& zoneName = iter.key();
            const label zonei = iter.val();

            cz[zonei] = new cellZone
            (
                zoneName,
                zoneCells[zonei].shrink(),
                zonei,
                mesh.cellZones()
            );
        }

        mesh.pointZones().resize(0);
        mesh.faceZones().resize(0);
        mesh.cellZones().resize(0);
        mesh.addZones(List<pointZone*>(), List<faceZone*>(), cz);
    }
}
