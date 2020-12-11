/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
    Copyright (C) 2019-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "sensitivityBezierFIIncompressible.H"
#include "addToRunTimeSelectionTable.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace incompressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(sensitivityBezierFI, 0);
addToRunTimeSelectionTable
(
    adjointSensitivity, sensitivityBezierFI, dictionary
);


// * * * * * * * * * * * Private  Member Functions  * * * * * * * * * * * * * //

void sensitivityBezierFI::read()
{
    // Laplace solution controls
    const dictionary dxdbDict = dict_.subOrEmptyDict("dxdbSolver");
    meshMovementIters_ = dxdbDict.getOrDefault<label>("iters", 1000);
    meshMovementResidualLimit_ =
        dxdbDict.getOrDefault<scalar>("tolerance", 1.e-07);

    // Read variables related to the adjoint eikonal solver
    FIBase::read();
}


tmp<volVectorField> sensitivityBezierFI::solveMeshMovementEqn
(
    const label iCP,
    const label idir
)
{
    read();
    tmp<volVectorField> tm(new volVectorField("m", dxdb_));
    volVectorField& m = tm.ref();

    // SOLVE FOR DXDB
    //~~~~~~~~~~~~~~~~
    // set boundary conditions
    for (const label patchI : sensitivityPatchIDs_)
    {
        // interpolate parameterization info to faces
        tmp<vectorField> tdxidXjFace = Bezier_.dxdbFace(patchI, iCP, idir);
        const vectorField& dxidXjFace = tdxidXjFace();

        m.boundaryFieldRef()[patchI] == dxidXjFace;
    }

    // iterate the adjoint to the eikonal equation
    for (label iter = 0; iter < meshMovementIters_; iter++)
    {
        Info<< "Mesh Movement Propagation(direction, CP), ("
            << idir << ", " << iCP << "), Iteration : "<< iter << endl;

        fvVectorMatrix mEqn
        (
            fvm::laplacian(m)
        );

        // Scalar residual = max(mEqn.solve().initialResidual());
        scalar residual = mag(mEqn.solve().initialResidual());

        Info<< "Max dxdb " << gMax(mag(m)()) << endl;

        mesh_.time().printExecutionTime(Info);

        // Check convergence
        if (residual < meshMovementResidualLimit_)
        {
            Info<< "\n***Reached dxdb convergence limit, iteration " << iter
                << "***\n\n";
            break;
        }
    }

    return tm;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sensitivityBezierFI::sensitivityBezierFI
(
    const fvMesh& mesh,
    const dictionary& dict,
    incompressibleVars& primalVars,
    incompressibleAdjointVars& adjointVars,
    objectiveManager& objectiveManager
)
:
    FIBase
    (
        mesh,
        dict,
        primalVars,
        adjointVars,
        objectiveManager
    ),
    //Bezier_(mesh, dict), // AJH Read locally?
    Bezier_(mesh, mesh.lookupObject<IOdictionary>("optimisationDict")),
    flowSens_(3*Bezier_.nBezier(), Zero),
    dSdbSens_(3*Bezier_.nBezier(), Zero),
    dndbSens_(3*Bezier_.nBezier(), Zero),
    dxdbDirectSens_(3*Bezier_.nBezier(), Zero),
    dVdbSens_(3*Bezier_.nBezier(), Zero),
    distanceSens_(3*Bezier_.nBezier(), Zero),
    optionsSens_(3*Bezier_.nBezier(), Zero),
    bcSens_(3*Bezier_.nBezier(), Zero),

    derivativesFolder_("optimisation"/type() + "Derivatives"),

    meshMovementIters_(-1),
    meshMovementResidualLimit_(1.e-7),
    dxdb_
    (
        variablesSet::autoCreateMeshMovementField
        (
            mesh,
            "mTilda",
            dimensionSet(dimLength)
        )
    )
{
    read();

    derivatives_ = scalarField(3*Bezier_.nBezier(), Zero),
    // Create folder to store sensitivities
    mkDir(derivativesFolder_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void sensitivityBezierFI::assembleSensitivities()
{
    // Adjoint to the eikonal equation
    autoPtr<volTensorField> distanceSensPtr(nullptr);
    if (includeDistance_)
    {
        // Solver equation
        eikonalSolver_->solve();

        // Allocate memory and compute grad(dxdb) multiplier
        distanceSensPtr.reset
        (
            createZeroFieldPtr<tensor>
            (
                mesh_,
                "distanceSensPtr",
                dimensionSet(0, 2, -3, 0, 0, 0, 0)
            )
        );
        distanceSensPtr() = eikonalSolver_->getFISensitivityTerm()().T();
    }

    const label nBezier = Bezier_.nBezier();
    const label nDVs = 3*nBezier;
    for (label iDV = 0; iDV < nDVs; iDV++)
    {
        label iCP = iDV%nBezier;
        label idir = iDV/nBezier;
        if
        (
            (idir == 0 && Bezier_.confineXmovement()[iCP])
         || (idir == 1 && Bezier_.confineYmovement()[iCP])
         || (idir == 2 && Bezier_.confineZmovement()[iCP])
        )
        {
            continue;
        }
        else
        {
            // Flow term
            // ~~~~~~~~~~~
            // compute dxdb and its grad
            tmp<volVectorField> tm = solveMeshMovementEqn(iCP, idir);
            const volVectorField& m = tm();
            volTensorField gradDxDb(fvc::grad(m, "grad(dxdb)"));

            flowSens_[iDV] =
                gSum
                (
                    (gradDxDb.primitiveField() && gradDxDbMult_.primitiveField())
                  * mesh_.V()
                );

            for (const label patchI : sensitivityPatchIDs_)
            {
                // Contribution from objective function
                // term from delta(n dS)/delta b and
                // term from delta(n)/delta b
                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                tmp<vectorField> tdSdb =
                    Bezier_.dndbBasedSensitivities(patchI, iCP, idir);
                const vectorField& dSdb = tdSdb();
                tmp<vectorField> tdndb =
                    Bezier_.dndbBasedSensitivities(patchI, iCP, idir, false);
                const vectorField& dndb = tdndb();
                dSdbSens_[iDV] += gSum(dSfdbMult_()[patchI] & dSdb);
                dndbSens_[iDV] += gSum(dnfdbMult_()[patchI] & dndb);

                // Contribution from objective function
                // term from delta( x ) / delta b
                // Only for objectives directly including
                // x, like moments
                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                tmp<vectorField> tdxdbFace =
                    Bezier_.dxdbFace(patchI, iCP, idir);
                const vectorField& dxdbFace = tdxdbFace();
                dxdbDirectSens_[iDV] +=
                    gSum(dxdbDirectMult_()[patchI] & dxdbFace);

                // Contribution from boundary conditions
                bcSens_[iDV] += gSum(bcDxDbMult_()[patchI] & dxdbFace);
            }

            // Contribution from delta (V) / delta b
            // For objectives defined as volume integrals only
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            dVdbSens_[iDV] =
                gSum
                (
                   divDxDbMult_
                 * fvc::div(m)().primitiveField()
                 * mesh_.V()
                );

            // Distance dependent term
            //~~~~~~~~~~~~~~~~~~~~~~~~~
            if (includeDistance_)
            {
                distanceSens_[iDV] =
                    gSum
                    (
                        (
                            distanceSensPtr().primitiveField()
                         && gradDxDb.primitiveField()
                        )
                       *mesh_.V()
                    );
            }

            // Terms from fvOptions
            optionsSens_[iDV] +=
                gSum((optionsDxDbMult_ & m.primitiveField())*mesh_.V());
        }

        // Sum contributions
        derivatives_ =
            flowSens_
          + dSdbSens_
          + dndbSens_
          + dxdbDirectSens_
          + dVdbSens_
          + distanceSens_
          + optionsSens_
          + bcSens_;
    }
}


void sensitivityBezierFI::clearSensitivities()
{
    flowSens_ = Zero;
    dSdbSens_ = Zero;
    dndbSens_ = Zero;
    dxdbDirectSens_ = Zero;
    dVdbSens_ = Zero;
    distanceSens_ = Zero;
    optionsSens_ = Zero;
    bcSens_ = Zero;

    FIBase::clearSensitivities();
}


void sensitivityBezierFI::write(const word& baseName)
{
    Info<< "Writing control point sensitivities to file" << endl;
    if (Pstream::master())
    {
        OFstream derivFile
        (
            derivativesFolder_/
                baseName + adjointVars_.solverName() + mesh_.time().timeName()
        );
        unsigned int widthDV = max(int(name(flowSens_.size()).size()), int(3));
        unsigned int width = IOstream::defaultPrecision() + 7;
        derivFile
            << setw(widthDV) << "#dv"        << " "
            << setw(width)   << "total"      << " "
            << setw(width)   << "flow"       << " "
            << setw(width)   << "dSdb"       << " "
            << setw(width)   << "dndb"       << " "
            << setw(width)   << "dxdbDirect" << " "
            << setw(width)   << "dVdb"       << " "
            << setw(width)   << "distance"   << " "
            << setw(width)   << "options"    << " "
            << setw(width)   << "dvdb"       << endl;
        const label nDVs = derivatives_.size();
        const label nBezier = Bezier_.nBezier();
        const boolListList& confineMovement = Bezier_.confineMovement();
        label lastActive(-1);

        for (label iDV = 0; iDV < nDVs; iDV++)
        {
            const label iCP(iDV%nBezier);
            const label idir(iDV/nBezier);
            if (!confineMovement[idir][iCP])
            {
                if (iDV!=lastActive + 1)
                {
                    derivFile << "\n";
                }
                lastActive = iDV;
                derivFile
                   << setw(widthDV) << iDV << " "
                   << setw(width) << derivatives_[iDV] << " "
                   << setw(width) << flowSens_[iDV] << " "
                   << setw(width) << dSdbSens_[iDV] << " "
                   << setw(width) << dndbSens_[iDV] << " "
                   << setw(width) << dxdbDirectSens_[iDV] << " "
                   << setw(width) << dVdbSens_[iDV] << " "
                   << setw(width) << distanceSens_[iDV] << " "
                   << setw(width) << optionsSens_[iDV] << " "
                   << setw(width) << bcSens_[iDV] << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
