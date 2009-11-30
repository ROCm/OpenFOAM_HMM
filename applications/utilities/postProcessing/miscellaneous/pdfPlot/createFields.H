    IOdictionary pdfDictionary
    (
        IOobject
        (
            "pdfDict",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    label nIntervals(readLabel(pdfDictionary.lookup("nIntervals")));

    label nSamples(readLabel(pdfDictionary.lookup("nSamples")));

    scalarField samples(nIntervals, 0);

