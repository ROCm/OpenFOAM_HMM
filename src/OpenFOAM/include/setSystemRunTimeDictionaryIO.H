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
    runTime,
    IOobject::MUST_READ_IF_MODIFIED,
    IOobject::NO_WRITE
);

if (dictPath.size())
{
    dictIO = IOobject
    (
        dictPath.expand(),
        runTime,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    );
}
