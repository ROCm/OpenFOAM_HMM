    argList::validArgs.clear();
    argList::noParallel();

    argList::validArgs.append("target root");
    argList::validArgs.append("target case");

    argList::validArgs.append("nLayers");
    argList::validArgs.append("overal thickness");


    argList::validOptions.insert("sourceRoot", "source root");
    argList::validOptions.insert("sourceCase", "source case");
    argList::validOptions.insert("sourcePatch", "source patch");

    argList::validOptions.insert("surface", "surface file");

    argList::validOptions.insert("mergeFaces", "");

    argList args(argc, argv);

    if (!args.check())
    {
         FatalError.exit();
    }


    fileName rootDirTarget(args.additionalArgs()[0]);
    fileName caseDirTarget(args.additionalArgs()[1]);

    label nLayers(readLabel(IStringStream(args.additionalArgs()[2])()));
    scalar thickness(readScalar(IStringStream(args.additionalArgs()[3])()));

    Info<< "Target: " << rootDirTarget << " " << caseDirTarget << endl;

    Info<< "Extruding layers:" << nl
        << "    number of layers " << nLayers << nl
        << "    overall thickness " << thickness << nl
        << endl;
