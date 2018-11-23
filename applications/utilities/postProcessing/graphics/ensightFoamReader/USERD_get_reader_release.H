int USERD_get_reader_release
(
    char release_number[Z_MAX_USERD_NAME]
)
{

#ifdef ENSIGHTDEBUG
    Info<< "Entering: USERD_get_reader_release" << endl;
#endif

    strncpy
    (
        release_number,
        #if OPENFOAM
        ("OpenFOAM-" + std::to_string(OPENFOAM)).c_str(),
        #else
        "OpenFOAM-unknown",
        #endif
        Z_MAX_USERD_NAME
    );

#ifdef ENSIGHTDEBUG
    Info<< "Leaving: USERD_get_reader_release" << endl;
#endif

    return Z_OK;
}
