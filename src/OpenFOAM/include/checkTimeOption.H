    if (args.optionFound("time"))
    {
        scalar timeValue = args.optionRead<scalar>("time");

        startTime = Time::findClosestTimeIndex(Times, timeValue);
        endTime = startTime + 1;
    }
