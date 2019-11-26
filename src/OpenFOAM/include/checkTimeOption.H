// Check -time and -latestTime options

if (args.found("time"))
{
    Foam::scalar timeValue = args.get<scalar>("time");

    startTime = Foam::Time::findClosestTimeIndex(Times, timeValue);
}

if (args.found("latestTime"))
{
    startTime = Times.size() - 1;
}
