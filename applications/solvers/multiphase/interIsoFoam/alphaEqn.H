    {
        // Note for AMR implementation:
        // At this point we have just entered the new time step,
        // the mesh has been refined and the alpha, phi and U contain
        // the field values at the beginning of this time step mapped
        // to the new mesh.

        // The above assumes that we are in firstIter() of the outer
        // corrector loop. If we are in any subsequent iter of this loop
        // the alpha1, U and phi will be overwritten with the new time step
        // values but still on the same mesh.


        if (pimple.firstIter())
        {
            // Note: This assumes moveMeshOuterCorrectors is false
            alpha10 = alpha1;
            U0 = U;
            phi0 = phi;
            advector.advect();

            #include "rhofs.H"
            rhoPhi = advector.getRhoPhi(rho1f, rho2f);
        }
        else
        {
            alpha1 = alpha10;
            // Temporarily setting U and phi to average of old and new value
            // Note: Illegal additions if mesh is topoChanging
            // (e.g. if moveMeshOuterCorrectors and AMR)
            U = 0.5*U0 + 0.5*U;
            phi = 0.5*phi0 + 0.5*phi;
            isoAdvection advector(alpha1, phi, U);
            advector.advect();
            #include "rhofs.H"
            rhoPhi = advector.getRhoPhi(rho1f, rho2f);
            // Resetting U and phi to the new value
            U = 2.0*U - U0;
            phi = 2.0*phi - phi0;
        }

        alpha2 = 1.0 - alpha1;
        mixture.correct();

    }
    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
