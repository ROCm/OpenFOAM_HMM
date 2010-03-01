    argList::validArgs.clear();

    argList::validArgs.append("master root");
    argList::validArgs.append("master case");
    argList::addOption("masterRegion", "name");

    argList::validArgs.append("root to add");
    argList::validArgs.append("case to add");
    argList::addOption("addRegion", "name");

    argList args(argc, argv);

    if (!args.check())
    {
         FatalError.exit();
    }

    fileName rootDirMaster = args[1];
    fileName caseDirMaster = args[2];
    word masterRegion = polyMesh::defaultRegion;
    args.optionReadIfPresent("masterRegion", masterRegion);

    fileName rootDirToAdd = args[3];
    fileName caseDirToAdd = args[4];
    word addRegion = polyMesh::defaultRegion;
    args.optionReadIfPresent("addRegion", addRegion);

    Info<< "Master:      " << rootDirMaster << " " << caseDirMaster
        << "  region " << masterRegion << nl
        << "mesh to add: " << rootDirToAdd << " " << caseDirToAdd
        << "  region " << addRegion << endl;

