/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
     \\/     M anipulation  | Copyright (C) 2017 IH-Cantabria
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

#include "multiphaseMangrovesTurbulenceModel.H"
#include "mathematicalConstants.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "fvmDdt.H"
#include "fvcGrad.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "surfaceInterpolate.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(multiphaseMangrovesTurbulenceModel, 0);
    addToRunTimeSelectionTable
    (
        option,
        multiphaseMangrovesTurbulenceModel,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::multiphaseMangrovesTurbulenceModel::multiphaseMangrovesTurbulenceModel
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    aZone_(),
    NZone_(),
    CkpZone_(),
    CepZone_(),
    CdZone_(),
    UName_("U"),
    kName_("k"),
    epsilonName_("epsilon")
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::multiphaseMangrovesTurbulenceModel::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);

    // Set alpha as zero-gradient
    const volScalarField& epsilon =
        mesh_.lookupObject<volScalarField>(epsilonName_);
    const volScalarField& k = mesh_.lookupObject<volScalarField>(kName_);

    volScalarField epsilonMangroves
    (
        IOobject
        (
            typeName + ":epsilonMangroves",
            mesh_.time().timeName(),
            mesh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("epsilonMangroves", dimMass/dimMass/dimTime, 0)
    );

    volScalarField kMangroves
    (
        IOobject
        (
            typeName + ":kMangroves",
            mesh_.time().timeName(),
            mesh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("kMangroves", dimTime/dimTime/dimTime, 0)
    );

    forAll(zoneIDs_, i)
    {
        const labelList& zones = zoneIDs_[i];

        forAll(zones, j)
        {
            const label zonei = zones[j];
            const cellZone& cz = mesh_.cellZones()[zonei];

            forAll(cz, kk)
            {
                const label celli = cz[kk];
                const scalar a = aZone_[i];
                const scalar N = NZone_[i];
                const scalar Ckp = CkpZone_[i];
                const scalar Cep = CepZone_[i];
                const scalar Cd = CdZone_[i];

                const vector& Uc = U[celli];

                const scalar f = Cd*a*N*mag(Uc);
                epsilonMangroves[celli] = Cep*f;
                kMangroves[celli] = Ckp*f;
            }
        }
    }

    if (eqn.psi().name() == epsilonName_)
    {
    	epsilonMangroves.correctBoundaryConditions();
	    fvMatrix<scalar> epsilonEqn
	    (
      	    - fvm::Sp(epsilonMangroves, epsilon)
    	);
	    eqn += epsilonEqn;
    }
    else if (eqn.psi().name() == kName_)
    {
    	kMangroves.correctBoundaryConditions();
	    fvMatrix<scalar> kEqn
	    (
      	    - fvm::Sp(kMangroves, k)
    	);
	    eqn += kEqn;
    }
}


bool Foam::fv::multiphaseMangrovesTurbulenceModel::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        if (coeffs_.found("epsilonNames"))
        {
            coeffs_.lookup("epsilonNames") >> fieldNames_;
        }
        else if (coeffs_.found("epsilon"))
        {
            word UName(coeffs_.lookup("epsilon"));
            fieldNames_ = wordList(1, UName);
        }
        else
        {
	        fieldNames_.setSize(2);
	        fieldNames_[0] = "epsilon";
    	    fieldNames_[1] = "k";
        }

        applied_.setSize(fieldNames_.size(), false);

        // Create the Mangroves models - 1 per region
        const dictionary& regionsDict(coeffs_.subDict("regions"));
        const wordList regionNames(regionsDict.toc());
        aZone_.setSize(regionNames.size(), 1);
        NZone_.setSize(regionNames.size(), 1);
        CkpZone_.setSize(regionNames.size(), 1);
        CepZone_.setSize(regionNames.size(), 1);
        CdZone_.setSize(regionNames.size(), 1);
        zoneIDs_.setSize(regionNames.size());

        forAll(zoneIDs_, i)
        {
            const word& regionName = regionNames[i];
            const dictionary& modelDict = regionsDict.subDict(regionName);

            const word& zoneName = modelDict.lookup("cellZone");
            zoneIDs_[i] = mesh_.cellZones().findIndices(zoneName);
            if (zoneIDs_[i].empty())
            {
                FatalErrorInFunction
                    << "Unable to find cellZone " << zoneName << nl
                    << "Valid cellZones are:" << mesh_.cellZones().names()
                    << exit(FatalError);
            }

            modelDict.lookup("a") >> aZone_[i];
            modelDict.lookup("N") >> NZone_[i];
            modelDict.lookup("Ckp") >> CkpZone_[i];
            modelDict.lookup("Cep") >> CepZone_[i];
            modelDict.lookup("Cd") >> CdZone_[i];
        }

        return true;
    }
    else
    {
        return false;
    }
}

void Foam::fv::multiphaseMangrovesTurbulenceModel::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);

    // Set alpha as zero-gradient
    const volScalarField& epsilon =
        mesh_.lookupObject<volScalarField>(epsilonName_);
    const volScalarField& k = mesh_.lookupObject<volScalarField>(kName_);

    volScalarField epsilonMangroves
    (
        IOobject
        (
            typeName + ":epsilonMangroves",
            mesh_.time().timeName(),
            mesh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("epsilonMangroves", dimMass/dimMass/dimTime, 0)
    );

    volScalarField kMangroves
    (
        IOobject
        (
            typeName + ":kMangroves",
            mesh_.time().timeName(),
            mesh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("kMangroves", dimTime/dimTime/dimTime, 0)
    );

    forAll(zoneIDs_, i)
    {
        const labelList& zones = zoneIDs_[i];

        forAll(zones, j)
        {
            const label zonei = zones[j];
            const cellZone& cz = mesh_.cellZones()[zonei];

            forAll(cz, kk)
            {
                const label celli = cz[kk];
                const scalar a = aZone_[i];
                const scalar N = NZone_[i];
                const scalar Ckp = CkpZone_[i];
                const scalar Cep = CepZone_[i];
                const scalar Cd = CdZone_[i];

                const vector& Uc = U[celli];

    epsilonMangroves[celli] = Cep * Cd * a * N * sqrt ( Uc.x()*Uc.x() + Uc.y()*Uc.y() + Uc.z()*Uc.z());
    kMangroves[celli] = Ckp * Cd * a * N * sqrt ( Uc.x()*Uc.x() + Uc.y()*Uc.y() + Uc.z()*Uc.z());
            }
        }
    }

    if (eqn.psi().name() == "epsilon")
    {
	Info<< " applying source to epsilon stuff.... " << endl;
    	epsilonMangroves.correctBoundaryConditions();
	fvMatrix<scalar> epsilonEqn
	(
      	    - fvm::Sp(rho*epsilonMangroves, epsilon)
    	);
	eqn += epsilonEqn;
    }
    else if (eqn.psi().name() == "k")
    {
	Info<< " applying source to k stuff.... " << endl;
    	kMangroves.correctBoundaryConditions();
	fvMatrix<scalar> kEqn
	(
      	    - fvm::Sp(rho*kMangroves, k)
    	);
	eqn += kEqn;
    }

}

void Foam::fv::multiphaseMangrovesTurbulenceModel::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    NotImplemented;
}

// ************************************************************************* //
