IOobject dictIO = IOobject::selectIO
(
    IOobject
    (
        dictName,
        runTime.system(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    ),
    args.opt<fileName>("dict", "")
);
