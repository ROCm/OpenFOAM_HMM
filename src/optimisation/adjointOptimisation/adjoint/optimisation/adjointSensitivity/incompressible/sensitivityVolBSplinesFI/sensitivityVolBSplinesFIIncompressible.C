/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "sensitivityVolBSplinesFIIncompressible.H"
#include "pointVolInterpolation.H"
#include "IOmanip.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace incompressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(sensitivityVolBSplinesFI, 0);
addToRunTimeSelectionTable
(
    adjointSensitivity,
    sensitivityVolBSplinesFI,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sensitivityVolBSplinesFI::sensitivityVolBSplinesFI
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
    volBSplinesBase_
    (
        const_cast<volBSplinesBase&>(volBSplinesBase::New(mesh))
    ),
    flowSens_(0),
    dSdbSens_(0),
    dndbSens_(0),
    dxdbDirectSens_(0),
    dVdbSens_(0),
    distanceSens_(0),
    optionsSens_(0),
    bcSens_(0),

    derivativesFolder_("optimisation"/type() + "Derivatives")
{
    // No boundary field pointers need to be allocated

    label nCPs = volBSplinesBase_.getTotalControlPointsNumber();
    derivatives_ = scalarField(3*nCPs, Zero);
    flowSens_ = vectorField(nCPs, Zero);
    dSdbSens_ = vectorField(nCPs, Zero);
    dndbSens_ = vectorField(nCPs, Zero);
    dxdbDirectSens_ = vectorField(nCPs, Zero);
    dVdbSens_ = vectorField(nCPs, Zero);
    distanceSens_ = vectorField(nCPs, Zero);
    optionsSens_ = vectorField(nCPs, Zero);
    bcSens_ = vectorField(nCPs, Zero);

    // Create folder to store sensitivities
    mkDir(derivativesFolder_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void sensitivityVolBSplinesFI::assembleSensitivities()
{
    read();

    // Interpolation engine
    pointVolInterpolation volPointInter(pointMesh::New(mesh_), mesh_);

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

    // Integration
    label passedCPs(0);
    PtrList<NURBS3DVolume>& boxes = volBSplinesBase_.boxesRef();
    forAll(boxes, iNURB)
    {
        const label nb(boxes[iNURB].getControlPoints().size());
        vectorField boxSensitivities(nb, Zero);

        vectorField dxdbSens = boxes[iNURB].computeControlPointSensitivities
        (
            dxdbDirectMult_(),
            sensitivityPatchIDs_.toc()
        );

        vectorField bcSens = boxes[iNURB].computeControlPointSensitivities
        (
            bcDxDbMult_(),
            sensitivityPatchIDs_.toc()
        );

        for (label cpI = 0; cpI < nb; cpI++)
        {
            label globalCP = passedCPs + cpI;

            // Parameterization info
            tmp<pointTensorField> dxdbI(boxes[iNURB].getDxDb(cpI));
            tmp<volTensorField> tvolDxDbI(volPointInter.interpolate(dxdbI));
            const volTensorField& volDxDbI = tvolDxDbI();

            // Chain rule used to get dx/db at cells
            // Gives practically the same results at a much higher CPU cost
            /*
            tmp<volTensorField> tvolDxDbI(boxes[iNURB].getDxCellsDb(cpI));
            volTensorField& volDxDbI = tvolDxDbI.ref();
            */

            // Gradient of parameterization info
            volVectorField temp
            (
                IOobject
                (
                    "dxdb",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                vector::zero
            );

            temp.replace(0, volDxDbI.component(0));
            temp.replace(1, volDxDbI.component(3));
            temp.replace(2, volDxDbI.component(6));
            volTensorField gradDxDb1(fvc::grad(temp));

            temp.replace(0, volDxDbI.component(1));
            temp.replace(1, volDxDbI.component(4));
            temp.replace(2, volDxDbI.component(7));
            volTensorField gradDxDb2(fvc::grad(temp));

            temp.replace(0, volDxDbI.component(2));
            temp.replace(1, volDxDbI.component(5));
            temp.replace(2, volDxDbI.component(8));
            volTensorField gradDxDb3(fvc::grad(temp));


            // Volume integral terms
            flowSens_[globalCP].x() = gSum
            (
                (gradDxDbMult_.primitiveField() && gradDxDb1.primitiveField())
               *mesh_.V()
            );
            flowSens_[globalCP].y() = gSum
            (
                (gradDxDbMult_.primitiveField() && gradDxDb2.primitiveField())
               *mesh_.V()
            );
            flowSens_[globalCP].z() = gSum
            (
                (gradDxDbMult_.primitiveField() && gradDxDb3.primitiveField())
               *mesh_.V()
            );

            // Contribution from objective function term from
            // delta( n dS ) / delta b and
            // delta ( x )   / delta b
            // for objectives directly depending on x
            for (const label patchI : sensitivityPatchIDs_)
            {
                tensorField dSdb
                (
                    boxes[iNURB].dndbBasedSensitivities(patchI, cpI)
                );
                dSdbSens_[globalCP] += gSum(dSfdbMult_()[patchI] & dSdb);
                tensorField dndb
                (
                    boxes[iNURB].dndbBasedSensitivities(patchI, cpI, false)
                );
                dndbSens_[globalCP] += gSum((dnfdbMult_()[patchI] & dndb));
            }

            // Contribution from delta (V) / delta b
            // For objectives defined as volume integrals only
            dVdbSens_[globalCP] +=
                gSum
                (
                    divDxDbMult_
                   *fvc::div(T(volDxDbI))().primitiveField()
                   *mesh_.V()
                );

            // Distance dependent term
            if (includeDistance_)
            {
                const tensorField& distSensInt =
                    distanceSensPtr().primitiveField();
                distanceSens_[globalCP].x() =
                    gSum
                    (
                        (distSensInt && gradDxDb1.primitiveField())*mesh_.V()
                    );
                distanceSens_[globalCP].y() =
                    gSum
                    (
                        (distSensInt && gradDxDb2.primitiveField())*mesh_.V()
                    );
                distanceSens_[globalCP].z() =
                    gSum
                    (
                        (distSensInt && gradDxDb3.primitiveField()) *mesh_.V()
                    );
            }

            // Terms from fvOptions
            optionsSens_[globalCP] +=
                gSum((optionsDxDbMult_ & volDxDbI.primitiveField())*mesh_.V());

            // dxdbSens storage
            dxdbDirectSens_[globalCP] = dxdbSens[cpI];

            // bcSens storage
            bcSens_[globalCP] = bcSens[cpI];

            boxSensitivities[cpI] =
                flowSens_[globalCP]
              + dSdbSens_[globalCP]
              + dndbSens_[globalCP]
              + dVdbSens_[globalCP]
              + distanceSens_[globalCP]
              + dxdbDirectSens_[globalCP]
              + optionsSens_[globalCP]
              + bcSens_[globalCP];
        }

        // Zero sensitivities in non-active design variables
        boxes[iNURB].boundControlPointMovement(boxSensitivities);

        // Transfer sensitivities to global list
        for (label cpI = 0; cpI < nb; cpI++)
        {
            label globalCP = passedCPs + cpI;
            derivatives_[3*globalCP] = boxSensitivities[cpI].x();
            derivatives_[3*globalCP + 1] = boxSensitivities[cpI].y();
            derivatives_[3*globalCP + 2] = boxSensitivities[cpI].z();
        }

        // Increment number of passed sensitivities
        passedCPs += nb;
    }

    // Zero non-active sensitivity components.
    // For consistent output only, does not affect optimisation
    volBSplinesBase_.boundControlPointMovement(flowSens_);
    volBSplinesBase_.boundControlPointMovement(dSdbSens_);
    volBSplinesBase_.boundControlPointMovement(dndbSens_);
    volBSplinesBase_.boundControlPointMovement(dVdbSens_);
    volBSplinesBase_.boundControlPointMovement(distanceSens_);
    volBSplinesBase_.boundControlPointMovement(dxdbDirectSens_);
    volBSplinesBase_.boundControlPointMovement(optionsSens_);
    volBSplinesBase_.boundControlPointMovement(bcSens_);
}


void sensitivityVolBSplinesFI::clearSensitivities()
{
    flowSens_ = vector::zero;
    dSdbSens_ = vector::zero;
    dndbSens_ = vector::zero;
    dxdbDirectSens_ = vector::zero;
    dVdbSens_ = vector::zero;
    distanceSens_ = vector::zero;
    optionsSens_ = vector::zero;
    bcSens_ = vector::zero;

    FIBase::clearSensitivities();
}


void sensitivityVolBSplinesFI::write(const word& baseName)
{
    Info<< "Writing control point sensitivities to file" << endl;
    if (Pstream::master())
    {
        OFstream derivFile
        (
             derivativesFolder_/
                baseName + adjointVars_.solverName() + mesh_.time().timeName()
        );
        unsigned int widthDV
        (
            max(int(Foam::name(flowSens_.size()).size()), int(3))
        );
        unsigned int width = IOstream::defaultPrecision() + 7;
        derivFile
            << setw(widthDV) << "#cp" << " "
            << setw(width) << "total::x" << " "
            << setw(width) << "total::y" << " "
            << setw(width) << "total::z" << " "
            << setw(width) << "flow::x" << " "
            << setw(width) << "flow::y" << " "
            << setw(width) << "flow::z" << " "
            << setw(width) << "dSdb::x" << " "
            << setw(width) << "dSdb::y" << " "
            << setw(width) << "dSdb::z" << " "
            << setw(width) << "dndb::x" << " "
            << setw(width) << "dndb::y" << " "
            << setw(width) << "dndb::z" << " "
            << setw(width) << "dxdbDirect::x" << " "
            << setw(width) << "dxdbDirect::y" << " "
            << setw(width) << "dxdbDirect::z" << " "
            << setw(width) << "dVdb::x" << " "
            << setw(width) << "dVdb::y" << " "
            << setw(width) << "dVdb::z" << " "
            << setw(width) << "distance::x" << " "
            << setw(width) << "distance::y" << " "
            << setw(width) << "distance::z" << " "
            << setw(width) << "options::x" << " "
            << setw(width) << "options::y" << " "
            << setw(width) << "options::z" << " "
            << setw(width) << "dvdb::x" << " "
            << setw(width) << "dvdb::y" << " "
            << setw(width) << "dvdb::z" << endl;

        label passedCPs(0);
        label lastActive(-1);
        PtrList<NURBS3DVolume>& boxes = volBSplinesBase_.boxesRef();
        forAll(boxes, iNURB)
        {
            label nb = boxes[iNURB].getControlPoints().size();
            const boolList& activeCPs = boxes[iNURB].getActiveCPs();
            for (label iCP = 0; iCP < nb; iCP++)
            {
                if (activeCPs[iCP])
                {
                    label globalCP = passedCPs + iCP;
                    if (globalCP!=lastActive + 1) derivFile << "\n";
                    lastActive = globalCP;

                    derivFile
                        << setw(widthDV) << globalCP << " "
                        << setw(width) << derivatives_[3*globalCP] << " "
                        << setw(width) << derivatives_[3*globalCP + 1] << " "
                        << setw(width) << derivatives_[3*globalCP + 2] << " "
                        << setw(width) << flowSens_[globalCP].x() << " "
                        << setw(width) << flowSens_[globalCP].y() << " "
                        << setw(width) << flowSens_[globalCP].z() << " "
                        << setw(width) << dSdbSens_[globalCP].x() << " "
                        << setw(width) << dSdbSens_[globalCP].y() << " "
                        << setw(width) << dSdbSens_[globalCP].z() << " "
                        << setw(width) << dndbSens_[globalCP].x() << " "
                        << setw(width) << dndbSens_[globalCP].y() << " "
                        << setw(width) << dndbSens_[globalCP].z() << " "
                        << setw(width) << dxdbDirectSens_[globalCP].x() << " "
                        << setw(width) << dxdbDirectSens_[globalCP].y() << " "
                        << setw(width) << dxdbDirectSens_[globalCP].z() << " "
                        << setw(width) << dVdbSens_[globalCP].x() << " "
                        << setw(width) << dVdbSens_[globalCP].y() << " "
                        << setw(width) << dVdbSens_[globalCP].z() << " "
                        << setw(width) << distanceSens_[globalCP].x() << " "
                        << setw(width) << distanceSens_[globalCP].y() << " "
                        << setw(width) << distanceSens_[globalCP].z() << " "
                        << setw(width) << optionsSens_[globalCP].x() << " "
                        << setw(width) << optionsSens_[globalCP].y() << " "
                        << setw(width) << optionsSens_[globalCP].z() << " "
                        << setw(width) << bcSens_[globalCP].x() << " "
                        << setw(width) << bcSens_[globalCP].y() << " "
                        << setw(width) << bcSens_[globalCP].z() << endl;
                }
            }
            passedCPs += nb;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
