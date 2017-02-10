fileName dictPath;
if (args.optionReadIfPresent("dict", dictPath))
{
    if (isDir(dictPath))
    {
        dictPath = dictPath / dictName;
    }
}

IOobject dictIO
(
    dictName,
    runTime.system(),
    mesh,
    IOobject::MUST_READ_IF_MODIFIED,
    IOobject::NO_WRITE
);

if (dictPath.size())
{
    dictIO = IOobject
    (
        dictPath.expand(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    );
}
