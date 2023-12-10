if (adjustTimeStep)
{
    runTime.setDeltaT(Foam::min(dtChem, maxDeltaT));
    Info<< "deltaT = " <<  runTime.deltaT().value() << endl;
}
