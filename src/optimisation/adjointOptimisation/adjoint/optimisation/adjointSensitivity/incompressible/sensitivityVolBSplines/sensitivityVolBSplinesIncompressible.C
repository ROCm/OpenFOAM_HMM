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

#include "sensitivityVolBSplinesIncompressible.H"
#include "addToRunTimeSelectionTable.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace incompressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(sensitivityVolBSplines, 0);
addToRunTimeSelectionTable
(
    adjointSensitivity,
    sensitivityVolBSplines,
    dictionary
);

// * * * * * * * * * * * Private  Member Functions  * * * * * * * * * * * * * //

void sensitivityVolBSplines::computeObjectiveContributions()
{
    if (includeObjective_)
    {
        label passedCPs = 0;
        PtrList<NURBS3DVolume>& boxes = volBSplinesBase_.boxesRef();
        forAll(boxes, iNURB)
        {
            label nb = boxes[iNURB].getControlPoints().size();
            for (label cpI = 0; cpI < nb; cpI++)
            {
                vector dSdbSensCP(Zero);
                vector dndbSensCP(Zero);
                for (const label patchI : sensitivityPatchIDs_)
                {
                    tensorField dSdb
                    (
                        boxes[iNURB].dndbBasedSensitivities(patchI, cpI)
                    );
                    dSdbSensCP += gSum(dSfdbMult_()[patchI] & dSdb);

                    tensorField dndb
                    (
                        boxes[iNURB].dndbBasedSensitivities
                        (
                            patchI,
                            cpI,
                            false
                        )
                    );
                    dndbSensCP += gSum((dnfdbMult_()[patchI] & dndb));
                }
                dSdbSens_[passedCPs + cpI] = dSdbSensCP;
                dndbSens_[passedCPs + cpI] = dndbSensCP;
            }
            passedCPs += nb;
        }
        volBSplinesBase_.boundControlPointMovement(dSdbSens_);
        volBSplinesBase_.boundControlPointMovement(dndbSens_);

        passedCPs = 0;
        forAll(boxes, iNURB)
        {
            vectorField sensDxDbDirect =
                boxes[iNURB].computeControlPointSensitivities
                (
                    dxdbDirectMult_(),
                    sensitivityPatchIDs_.toc()
                );

            // Transfer to global list
            forAll(sensDxDbDirect, cpI)
            {
                dxdbDirectSens_[passedCPs + cpI] = sensDxDbDirect[cpI];
            }
            passedCPs += sensDxDbDirect.size();
        }
        volBSplinesBase_.boundControlPointMovement(dxdbDirectSens_);
    }
}


void sensitivityVolBSplines::computeBCContributions()
{
    label passedCPs = 0;
    PtrList<NURBS3DVolume>& boxes = volBSplinesBase_.boxesRef();
    forAll(boxes, iNURB)
    {
        vectorField sensBcsDxDb =
            boxes[iNURB].computeControlPointSensitivities
            (
                bcDxDbMult_(),
                sensitivityPatchIDs_.toc()
            );

        // Transfer to global list
        forAll(sensBcsDxDb, cpI)
        {
            bcSens_[passedCPs + cpI] = sensBcsDxDb[cpI];
        }
        passedCPs += sensBcsDxDb.size();
    }
    volBSplinesBase_.boundControlPointMovement(bcSens_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sensitivityVolBSplines::sensitivityVolBSplines
(
    const fvMesh& mesh,
    const dictionary& dict,
    incompressibleVars& primalVars,
    incompressibleAdjointVars& adjointVars,
    objectiveManager& objectiveManager
)
:
    SIBase
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
    bcSens_(0),

    derivativesFolder_("optimisation"/type() + "Derivatives")
{
    // No boundary field pointers need to be allocated
    const label nCPs(volBSplinesBase_.getTotalControlPointsNumber());
    derivatives_ = scalarField(3*nCPs, Zero);
    flowSens_ = vectorField(nCPs, Zero);
    dSdbSens_ = vectorField(nCPs, Zero);
    dndbSens_ = vectorField(nCPs, Zero);
    dxdbDirectSens_ = vectorField(nCPs, Zero);
    bcSens_ = vectorField(nCPs, Zero);

    // Create folder to store sensitivities
    mkDir(derivativesFolder_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void sensitivityVolBSplines::assembleSensitivities()
{
    // Assemble the sensitivity map
    // Solves for the post-processing equations and adds their contribution to
    // the sensitivity map
    surfaceSensitivity_.assembleSensitivities();

    // Finalise sensitivities including dxFace/db
    const boundaryVectorField& faceSens =
        surfaceSensitivity_.getWallFaceSensVecBoundary();

    label passedCPs(0);
    PtrList<NURBS3DVolume>& boxes = volBSplinesBase_.boxesRef();
    forAll(boxes, iNURB)
    {
        vectorField sens =
            boxes[iNURB].computeControlPointSensitivities
            (
                faceSens,
                sensitivityPatchIDs_.toc()
            );
        // Transfer to global list
        forAll(sens, cpI)
        {
            flowSens_[passedCPs + cpI] = sens[cpI];
        }
        passedCPs += sens.size();
    }
    volBSplinesBase_.boundControlPointMovement(flowSens_);

    // Contribution from objective function
    // Note:
    // includeObjectiveContribution has to be set to false (false by default)
    // in surfaceSensitivity, in order to avoid computing this term twice.
    // Optionally avoided altogether if includeObjectiveContribution is set to
    // false for sensitivityVolBSplines
    computeObjectiveContributions();

    computeBCContributions();

    // Transform sensitivites to scalarField in order to cooperate with
    // updateMethod
    forAll(flowSens_, cpI)
    {
        derivatives_[3*cpI] =
            flowSens_[cpI].x()
          + dSdbSens_[cpI].x()
          + dndbSens_[cpI].x()
          + dxdbDirectSens_[cpI].x()
          + bcSens_[cpI].x();
        derivatives_[3*cpI + 1] =
            flowSens_[cpI].y()
          + dSdbSens_[cpI].y()
          + dndbSens_[cpI].y()
          + dxdbDirectSens_[cpI].y()
          + bcSens_[cpI].y();
        derivatives_[3*cpI + 2] =
            flowSens_[cpI].z()
          + dSdbSens_[cpI].z()
          + dndbSens_[cpI].z()
          + dxdbDirectSens_[cpI].z()
          + bcSens_[cpI].z();
    }
}


void sensitivityVolBSplines::clearSensitivities()
{
    flowSens_ = vector::zero;
    dSdbSens_ = vector::zero;
    dndbSens_ = vector::zero;
    dxdbDirectSens_ = vector::zero;
    bcSens_ = vector::zero;

    SIBase::clearSensitivities();
}


void sensitivityVolBSplines::write(const word& baseName)
{
    Info<< "Writing control point sensitivities to file" << endl;
    // Write sensitivity map
    SIBase::write(baseName);
    // Write control point sensitivities
    if (Pstream::master())
    {
        OFstream derivFile
        (
             derivativesFolder_/
                 baseName + adjointVars_.solverName() + mesh_.time().timeName()
        );
        unsigned int widthDV =
            max(int(Foam::name(derivatives_.size()).size()), int(3));
        unsigned int width = IOstream::defaultPrecision() + 7;
        derivFile
            << setw(widthDV) << "#cp" << " "
            << setw(width) << "total::x"<< " "
            << setw(width) << "total::y"<< " "
            << setw(width) << "total::z"<< " "
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
                    if (globalCP!=lastActive + 1)
                    {
                        derivFile << "\n";
                    }
                    lastActive = globalCP;

                    derivFile
                       << setw(widthDV) << globalCP << " "
                       << setw(width) << derivatives_[3*globalCP]  << " "
                       << setw(width) << derivatives_[3*globalCP + 1]  << " "
                       << setw(width) << derivatives_[3*globalCP + 2]  << " "
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
                       << setw(width) << bcSens_[globalCP].x() << " "
                       << setw(width) << bcSens_[globalCP].y() << " "
                       << setw(width) << bcSens_[globalCP].z()
                       << endl;
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
