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

const word setFormat(propsDict.getOrDefault<word>("setFormat", "vtk"));
