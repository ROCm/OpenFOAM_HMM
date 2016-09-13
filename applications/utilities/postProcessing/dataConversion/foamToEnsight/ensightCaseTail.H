if (Pstream::master())
{
    ensightCaseFile.setf(ios_base::scientific, ios_base::floatfield);
    ensightCaseFile.precision(5);

    ensightCaseFile << nl << "TIME" << nl
        << "time set:                      " << 1 << nl
        << "number of steps:               " << nTimeSteps << nl
        << "filename start number:         " << 0 << nl
        << "filename increment:            " << 1 << nl;

    ensightCaseFile << "time values:" << nl;

    label count = 0;
    scalar Tcorr = 0.0;
    if (timeDirs[0].value() < 0)
    {
        Tcorr = -timeDirs[0].value();
        Info<< "Correcting time values. Adding " << Tcorr << endl;
    }

    forAll(timeDirs, n)
    {
        ensightCaseFile << setw(12) << timeDirs[n].value() + Tcorr << " ";

        if (++count % 6 == 0)
        {
            ensightCaseFile << nl;
        }
    }

    ensightCaseFile << nl;
}
