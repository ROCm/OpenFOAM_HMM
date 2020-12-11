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

#include "sensitivityBezierIncompressible.H"
#include "addToRunTimeSelectionTable.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace incompressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(sensitivityBezier, 0);
addToRunTimeSelectionTable
(
    adjointSensitivity,
    sensitivityBezier,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sensitivityBezier::sensitivityBezier
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
    //Bezier_(mesh, dict), // AJH Read locally?
    Bezier_(mesh, mesh.lookupObject<IOdictionary>("optimisationDict")),
    sens_(Bezier_.nBezier(), Zero),
    flowSens_(Bezier_.nBezier(), Zero),
    dSdbSens_(Bezier_.nBezier(), Zero),
    dndbSens_(Bezier_.nBezier(), Zero),
    dxdbDirectSens_(Bezier_.nBezier(), Zero),
    bcSens_(Bezier_.nBezier(), Zero),
    derivativesFolder_("optimisation"/type() + "Derivatives")
{
    derivatives_ = scalarField(3*Bezier_.nBezier(), Zero);
    // Create folder to store sensitivities
    mkDir(derivativesFolder_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void sensitivityBezier::assembleSensitivities()
{
    // Assemble the sensitivity map
    // Solves for the post-processing equations and adds their contribution to
    // the sensitivity map
    surfaceSensitivity_.assembleSensitivities();

    forAll(sens_, iCP)
    {
        // Face-based summation. More robust since the result is independent of
        // the number of processors (does not hold for a point-based summation)
        for (const label patchI : sensitivityPatchIDs_)
        {
            // Interpolate parameterization info to faces
            tmp<tensorField> tdxidXj = Bezier_.dxdbFace(patchI, iCP);
            const tensorField& dxidXj = tdxidXj();

            // Patch sensitivity map
            const vectorField& patchSensMap =
                surfaceSensitivity_.getWallFaceSensVecBoundary()[patchI];
            flowSens_[iCP] += gSum(patchSensMap & dxidXj);

            if (includeObjective_)
            {
                // Contribution from objective function
                // term from delta( n dS ) / delta b and
                // term from delta( n    ) / delta b
                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                tmp<tensorField> tdSdb
                (
                    Bezier_.dndbBasedSensitivities(patchI, iCP)
                );
                const tensorField& dSdb = tdSdb();
                dSdbSens_[iCP] += gSum(dSfdbMult_()[patchI] & dSdb);

                tmp<tensorField> tdndb
                (
                    Bezier_.dndbBasedSensitivities(patchI, iCP, false)
                );
                const tensorField& dndb = tdndb();
                dndbSens_[iCP] += gSum((dnfdbMult_()[patchI] & dndb));

                // Contribution from objective function
                // term from delta( x ) / delta b
                // Only for objectives directly including
                // x, like moments
                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                dxdbDirectSens_[iCP] +=
                    gSum((dxdbDirectMult_()[patchI] & dxidXj));
            }

            // Sensitivities from boundary conditions
            bcSens_[iCP] += gSum(bcDxDbMult_()[patchI] & dxidXj);
        }
    }
    sens_ = flowSens_ + dSdbSens_ + dndbSens_ + dxdbDirectSens_ + bcSens_;

    // Transform sensitivities to scalarField in order to cooperate with
    // updateMethod
    label nBezier = Bezier_.nBezier();
    forAll(sens_, cpI)
    {
        derivatives_[cpI] = sens_[cpI].x();
        derivatives_[cpI + nBezier] = sens_[cpI].y();
        derivatives_[cpI + 2*nBezier] = sens_[cpI].z();
        const boolList& confineXmovement = Bezier_.confineXmovement();
        const boolList& confineYmovement = Bezier_.confineYmovement();
        const boolList& confineZmovement = Bezier_.confineZmovement();
        if (confineXmovement[cpI])
        {
            derivatives_[cpI] *= scalar(0);
            flowSens_[cpI].x() = Zero;
            dSdbSens_[cpI].x() = Zero;
            dndbSens_[cpI].x() = Zero;
            dxdbDirectSens_[cpI].x() = Zero;
            bcSens_[cpI].x() = Zero;
        }
        if (confineYmovement[cpI])
        {
            derivatives_[cpI + nBezier] *= scalar(0);
            flowSens_[cpI].y() = Zero;
            dSdbSens_[cpI].y() = Zero;
            dndbSens_[cpI].y() = Zero;
            dxdbDirectSens_[cpI].y() = Zero;
            bcSens_[cpI].y() = Zero;
        }
        if (confineZmovement[cpI])
        {
            derivatives_[cpI + 2*nBezier] *= scalar(0);
            flowSens_[cpI].z() = Zero;
            dSdbSens_[cpI].z() = Zero;
            dndbSens_[cpI].z() = Zero;
            dxdbDirectSens_[cpI].z() = Zero;
            bcSens_[cpI].z() = Zero;
        }
    }
}


void sensitivityBezier::clearSensitivities()
{
    sens_ = Zero;
    flowSens_ = Zero;
    dSdbSens_ = Zero;
    dndbSens_ = Zero;
    dxdbDirectSens_ = Zero;
    bcSens_ = Zero;

    SIBase::clearSensitivities();
}


void sensitivityBezier::write(const word& baseName)
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
        unsigned int widthDV = max(int(name(sens_.size()).size()), int(3));
        unsigned int width = IOstream::defaultPrecision() + 7;
        derivFile
            << setw(widthDV) << "#dv" << " "
            << setw(width) << "total" << " "
            << setw(width) << "flow" << " "
            << setw(width) << "dSdb" << " "
            << setw(width) << "dndb" << " "
            << setw(width) << "dxdbDirect" << " "
            << setw(width) << "dvdb" << endl;
        label nDV = derivatives_.size();
        label nBezier = Bezier_.nBezier();
        const boolListList& confineMovement = Bezier_.confineMovement();
        label lastActive(-1);
        for (label iDV = 0; iDV < nDV; iDV++)
        {
            label iCP = iDV%nBezier;
            label idir = iDV/nBezier;
            if (!confineMovement[idir][iCP])
            {
                if (iDV!=lastActive + 1) derivFile << "\n";
                lastActive = iDV;
                derivFile
                    << setw(widthDV) << iDV << " "
                    << setw(width) << derivatives_[iDV] << " "
                    << setw(width) << flowSens_[iCP].component(idir) << " "
                    << setw(width) << dSdbSens_[iCP].component(idir) << " "
                    << setw(width) << dndbSens_[iCP].component(idir) << " "
                    << setw(width) << dxdbDirectSens_[iCP].component(idir) << " "
                    << setw(width) << bcSens_[iCP].component(idir)
                    << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
