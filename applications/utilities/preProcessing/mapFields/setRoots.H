    argList::validArgs.clear();

    argList::validOptions.insert("source", "dir");

    argList::validOptions.erase(argList::validOptions.find("parallel"));
    argList::validOptions.insert("parallelSource", "");
    argList::validOptions.insert("parallelTarget", "");
    argList::validOptions.insert("consistent", "");

    argList args(argc, argv);

    if (!args.check())
    {
         FatalError.exit();
    }

    fileName rootDirTarget(args.rootPath());
    fileName caseDirTarget(args.globalCaseName());

    if (!args.options().found("source"))
    {
        FatalErrorIn(args.executable())
            << "Essential argument [-source dir] not provided"
            << exit(FatalError);
    }

    fileName casePath(args.options()["source"]);
    fileName rootDirSource = casePath.path();
    fileName caseDirSource = casePath.name();

    Info<< "Source: " << rootDirSource << " " << caseDirSource << nl
        << "Target: " << rootDirTarget << " " << caseDirTarget << endl;

    bool parallelSource = false;
    if (args.options().found("parallelSource"))
    {
        parallelSource = true;
    }

    bool parallelTarget = false;
    if (args.options().found("parallelTarget"))
    {
        parallelTarget = true;
    }

    bool consistent = false;
    if (args.options().found("consistent"))
    {
        consistent = true;
    }
