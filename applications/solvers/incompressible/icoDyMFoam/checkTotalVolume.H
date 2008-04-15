    scalar newTotalVolume = sum(mesh.V());
    scalar totalVolRatio = newTotalVolume/totalVolume;

    Info << "Total volume change: " << totalVolRatio - 1 << endl;

    totalVolume = newTotalVolume;
