{
    // Semi-implicit mass transfer for species
    // Initilize dmdt for alpha Eq's for mass transfers driven by species
    autoPtr<phaseSystem::massTransferTable>
        massTransferPtr(fluid.massTransfer(T));

    //phaseSystem::massTransferTable& massTransfer(massTransferPtr());

    forAllIter(UPtrList<phaseModel>, fluid.phases(), iter)
    {
        phaseModel& phase = iter();
        PtrList<volScalarField>& Y = phase.Y();
        //const surfaceScalarField& alphaPhi = phase.alphaPhi();

        if (!Y.empty())
        {
            const volScalarField& alpha = phase;

            label inertIndex = -1;
            volScalarField Yt(0.0*Y[0]);

            forAll(Y, i)
            {
                tmp<fvScalarMatrix> YiEqn(phase.YiEqn(Y[i]));

                if (YiEqn.valid())
                {
                    YiEqn.ref() =
                    (
                        YiEqn()
                      - fvm::laplacian
                        (
                            alpha*turbulence->nuEff(),
                            Y[i]
                        )
                    ==
                    //  (*massTransfer[Y[i].name()])(/phase.rho()
                        fvc::ddt(alpha)
                       *pos
                        (
                           fluid.dmdtYi(Y[i].name())
                        - dimensionedScalar("zero", dimDensity/dimTime, 1e-3)
                        )
                        //fluid.dmdtYi(Y[i].name())/phase.rho()
                       //explicit mass sources (P or T)
                    );

                    YiEqn->relax();
                    YiEqn->solve(mesh.solver("Yi"));
                    Y[i].max(0.0);
                    Y[i].min(1.0);
                    Yt += Y[i];
                }
                else
                {
                    inertIndex = i;
                }

//                 if (mesh.time().outputTime())
//                 {
//                     volScalarField dmdtYi("dmdtYi", pos(fluid.dmdtYi(Y[i].name())));
//                     dmdtYi.write();
//                 }

                Info << "Min/Max : " << min(Y[i]) << " " << max(Y[i]) << endl;
                Info<< "Max dmdtYi : "
                    << max(fluid.dmdtYi(Y[i].name())().internalField()) << endl;
                Info<< "Min dmdtYi : "
                    << min(fluid.dmdtYi(Y[i].name())().internalField()) << endl;
            }
            Y[inertIndex] = scalar(1) - Yt;
            Y[inertIndex].max(0.0);
        }
    }
}

