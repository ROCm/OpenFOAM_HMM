    // Momentum predictor

    #ifdef USE_ROCTX
    roctxRangePush("UEqn");
    #endif

    #ifdef USE_ROCTX
    roctxRangePush("UEqn-MRF.correctBoundaryVelocity");
    #endif

    MRF.correctBoundaryVelocity(U);

    #ifdef USE_ROCTX
    roctxRangePop();
    #endif

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
    );

    fvVectorMatrix& UEqn = tUEqn.ref();

    #ifdef USE_ROCTX
    roctxRangePush("UEqn-UEqn.relax");
    #endif

    UEqn.relax();

    #ifdef USE_ROCTX
    roctxRangePop();
    #endif

    #ifdef USE_ROCTX
    roctxRangePush("UEqn-fvOptions.constrain");
    #endif

    fvOptions.constrain(UEqn);

    #ifdef USE_ROCTX
    roctxRangePop();
    #endif


    if (simple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvOptions.correct(U);
    }

    #ifdef USE_ROCTX
    roctxRangePop();
    #endif
