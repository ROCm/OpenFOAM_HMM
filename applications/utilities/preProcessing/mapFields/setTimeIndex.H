    label sourceTimeIndex = runTimeSource.timeIndex();
    instantList sourceTimes = runTimeSource.times();
    if (args.options().found("sourceTime"))
    {
        if ((args.options()["sourceTime"]) == "latestTime")
        {
            sourceTimeIndex = sourceTimes.size() - 1;
        }
        else
        {
            scalar sourceTime =
                readScalar(IStringStream(args.options()["sourceTime"])());
            sourceTimeIndex = getTimeIndex(sourceTimes, sourceTime);
        }
    }
