IOdictionary propsDict
(
    IOobject
    (
        "particleTrackProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED
    )
);

const word cloudName(propsDict.get<word>("cloud"));

const label sampleFrequency(propsDict.get<label>("sampleFrequency"));

const label maxPositions(propsDict.get<label>("maxPositions"));

const label maxTracks(propsDict.getOrDefault<label>("maxTracks", -1));

const word setFormat(propsDict.getOrDefault<word>("setFormat", "vtk"));

const wordRes fieldNames(propsDict.getOrDefault<wordRes>("fields", wordRes()));

const word UName(propsDict.getOrDefault<word>("U", "U"));

const dictionary formatOptions = propsDict.subOrEmptyDict("formatOptions");
