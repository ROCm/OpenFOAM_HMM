    #include "readTimeControls.H"
    #include "alphaControls.H"

    Switch implicitPhasePressure
    (
        alphaControls.lookupOrDefault<Switch>("implicitPhasePressure", false)
    );
