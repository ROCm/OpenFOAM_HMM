    instantList sourceTimes = runTimeSource.times();
    label sourceTimeIndex = runTimeSource.timeIndex();
    if (args.optionFound("sourceTime"))
    {
        if (args.option("sourceTime") == "latestTime")
        {
            sourceTimeIndex = sourceTimes.size() - 1;
        }
        else
        {
            scalar sourceTime = args.optionRead<scalar>("sourceTime");
            sourceTimeIndex = getTimeIndex(sourceTimes, sourceTime);
        }
    }
    else
    {
        scalar sourceTime = runTimeTarget.time().value();
        sourceTimeIndex = getTimeIndex(sourceTimes, sourceTime);
    }
