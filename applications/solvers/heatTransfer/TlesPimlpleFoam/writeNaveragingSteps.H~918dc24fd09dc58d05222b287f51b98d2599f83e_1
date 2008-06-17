    if (runTime.outputTime())
    {
        OFstream nAveragingStepsFile
        (
            runTime.path()/runTime.timeName()/"uniform"/"nAveragingSteps.raw"
        );

        if (nAveragingStepsFile.good())
        {
            nAveragingStepsFile << nAveragingSteps << endl;
        }
        else
        {
            FatalErrorIn(args.executable())
                << "Cannot open file " << nAveragingStepsFile.name()
                << abort(FatalError);
        }
    }
