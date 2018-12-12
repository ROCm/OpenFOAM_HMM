IOobject dictIO = IOobject::selectIO
(
    IOobject
    (
        dictName,
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    ),
    args.opt<fileName>("dict", "")
);
