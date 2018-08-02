IOobject dictIO = IOobject::selectIO
(
    IOobject
    (
        dictName,
        runTime.system(),
        runTime,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    ),
    args.lookupOrDefault<fileName>("dict", "")
);
