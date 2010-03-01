    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("source dir");

    argList::addOption("sourceTime", "scalar");
    argList::addBoolOption("parallelSource");
    argList::addBoolOption("parallelTarget");
    argList::addBoolOption("consistent");

    argList args(argc, argv);

    if (!args.check())
    {
         FatalError.exit();
    }

    fileName rootDirTarget(args.rootPath());
    fileName caseDirTarget(args.globalCaseName());

    const fileName casePath = args[1];
    const fileName rootDirSource = casePath.path();
    const fileName caseDirSource = casePath.name();

    Info<< "Source: " << rootDirSource << " " << caseDirSource << nl
        << "Target: " << rootDirTarget << " " << caseDirTarget << endl;

    const bool parallelSource = args.optionFound("parallelSource");
    const bool parallelTarget = args.optionFound("parallelTarget");
    const bool consistent     = args.optionFound("consistent");

