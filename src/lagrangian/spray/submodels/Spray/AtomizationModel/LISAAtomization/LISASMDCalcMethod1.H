{
    // calculate the new diameter with a Rosin Rammler distribution

    scalar minValue = min(d, dD/10.0);
    scalar maxValue = dD;

    if(maxValue - minValue < SMALL)
    {
        minValue = d/10.0;
    }

    scalar range = maxValue - minValue;

    scalar y = 0;

    bool success = false;

    while(!success)
    {

        x = minValue + range*rndGen.scalar01();
        y = rndGen.scalar01();
        scalar p = 0.0;
        scalar nExp = 1;
        scalar xx = pow(x/dD, nExp);

        p = xx*exp(-xx);
        if (y<p)
        {
            success = true;
        }
    }
}
