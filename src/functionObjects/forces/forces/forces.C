/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "forces.H"
#include "fvcGrad.H"
#include "porosityModel.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "cartesianCS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(forces, 0);
    addToRunTimeSelectionTable(functionObject, forces, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::forces::setCoordinateSystem
(
    const dictionary& dict,
    const word& e3Name,
    const word& e1Name
)
{
    coordSysPtr_.clear();

    point origin(Zero);
    if (dict.readIfPresent<point>("CofR", origin))
    {
        const vector e3 = e3Name == word::null ?
            vector(0, 0, 1) : dict.get<vector>(e3Name);
        const vector e1 = e1Name == word::null ?
            vector(1, 0, 0) : dict.get<vector>(e1Name);

        coordSysPtr_.reset(new coordSystem::cartesian(origin, e3, e1));
    }
    else
    {
        // The 'coordinateSystem' sub-dictionary is optional,
        // but enforce use of a cartesian system if not found.

        if (dict.found(coordinateSystem::typeName_()))
        {
            // New() for access to indirect (global) coordinate system
            coordSysPtr_ =
                coordinateSystem::New
                (
                    obr_,
                    dict,
                    coordinateSystem::typeName_()
                );
        }
        else
        {
            coordSysPtr_.reset(new coordSystem::cartesian(dict));
        }
    }
}


Foam::volVectorField& Foam::functionObjects::forces::force()
{
    auto* forcePtr = mesh_.getObjectPtr<volVectorField>(scopedName("force"));

    if (!forcePtr)
    {
        forcePtr = new volVectorField
        (
            IOobject
            (
                scopedName("force"),
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector(dimForce, Zero)
        );

        mesh_.objectRegistry::store(forcePtr);
    }

    return *forcePtr;
}


Foam::volVectorField& Foam::functionObjects::forces::moment()
{
    auto* momentPtr = mesh_.getObjectPtr<volVectorField>(scopedName("moment"));

    if (!momentPtr)
    {
        momentPtr = new volVectorField
        (
            IOobject
            (
                scopedName("moment"),
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector(dimForce*dimLength, Zero)
        );

        mesh_.objectRegistry::store(momentPtr);
    }

    return *momentPtr;
}


void Foam::functionObjects::forces::initialise()
{
    if (initialised_)
    {
        return;
    }

    if (directForceDensity_)
    {
        if (!foundObject<volVectorField>(fDName_))
        {
            FatalErrorInFunction
                << "Could not find " << fDName_ << " in database"
                << exit(FatalError);
        }
    }
    else
    {
        if
        (
            !foundObject<volVectorField>(UName_)
         || !foundObject<volScalarField>(pName_)
        )
        {
            FatalErrorInFunction
                << "Could not find U: " << UName_
                << " or p:" << pName_ << " in database"
                << exit(FatalError);
        }

        if (rhoName_ != "rhoInf" && !foundObject<volScalarField>(rhoName_))
        {
            FatalErrorInFunction
                << "Could not find rho:" << rhoName_ << " in database"
                << exit(FatalError);
        }
    }

    initialised_ = true;
}


void Foam::functionObjects::forces::reset()
{
    sumPatchForcesP_ = Zero;
    sumPatchForcesV_ = Zero;
    sumPatchMomentsP_ = Zero;
    sumPatchMomentsV_ = Zero;

    sumInternalForces_ = Zero;
    sumInternalMoments_ = Zero;

    auto& force = this->force();
    auto& moment = this->moment();
    force == dimensionedVector(force.dimensions(), Zero);
    moment == dimensionedVector(moment.dimensions(), Zero);
}


Foam::tmp<Foam::volSymmTensorField>
Foam::functionObjects::forces::devRhoReff() const
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    if (foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        const auto& turb =
            lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);

        return turb.devRhoReff();
    }
    else if (foundObject<icoTurbModel>(icoTurbModel::propertiesName))
    {
        const auto& turb =
            lookupObject<icoTurbModel>(icoTurbModel::propertiesName);

        return rho()*turb.devReff();
    }
    else if (foundObject<fluidThermo>(fluidThermo::dictName))
    {
        const auto& thermo = lookupObject<fluidThermo>(fluidThermo::dictName);

        const auto& U = lookupObject<volVectorField>(UName_);

        return -thermo.mu()*dev(twoSymm(fvc::grad(U)));
    }
    else if (foundObject<transportModel>("transportProperties"))
    {
        const auto& laminarT =
            lookupObject<transportModel>("transportProperties");

        const auto& U = lookupObject<volVectorField>(UName_);

        return -rho()*laminarT.nu()*dev(twoSymm(fvc::grad(U)));
    }
    else if (foundObject<dictionary>("transportProperties"))
    {
        const auto& transportProperties =
            lookupObject<dictionary>("transportProperties");

        const dimensionedScalar nu("nu", dimViscosity, transportProperties);

        const auto& U = lookupObject<volVectorField>(UName_);

        return -rho()*nu*dev(twoSymm(fvc::grad(U)));
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for viscous stress calculation"
            << exit(FatalError);

        return volSymmTensorField::null();
    }
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::forces::mu() const
{
    if (foundObject<fluidThermo>(basicThermo::dictName))
    {
        const auto& thermo = lookupObject<fluidThermo>(basicThermo::dictName);

        return thermo.mu();
    }
    else if (foundObject<transportModel>("transportProperties"))
    {
        const auto& laminarT =
            lookupObject<transportModel>("transportProperties");

        return rho()*laminarT.nu();
    }
    else if (foundObject<dictionary>("transportProperties"))
    {
        const auto& transportProperties =
             lookupObject<dictionary>("transportProperties");

        const dimensionedScalar nu("nu", dimViscosity, transportProperties);

        return rho()*nu;
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for dynamic viscosity calculation"
            << exit(FatalError);

        return volScalarField::null();
    }
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::forces::rho() const
{
    if (rhoName_ == "rhoInf")
    {
        return tmp<volScalarField>::New
        (
            IOobject
            (
                "rho",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar(dimDensity, rhoRef_)
        );
    }

    return (lookupObject<volScalarField>(rhoName_));
}


Foam::scalar Foam::functionObjects::forces::rho(const volScalarField& p) const
{
    if (p.dimensions() == dimPressure)
    {
        return 1;
    }

    if (rhoName_ != "rhoInf")
    {
        FatalErrorInFunction
            << "Dynamic pressure is expected but kinematic is provided."
            << exit(FatalError);
    }

    return rhoRef_;
}


void Foam::functionObjects::forces::addToPatchFields
(
    const label patchi,
    const vectorField& Md,
    const vectorField& fP,
    const vectorField& fV
)
{
    sumPatchForcesP_ += sum(fP);
    sumPatchForcesV_ += sum(fV);
    force().boundaryFieldRef()[patchi] += fP + fV;

    const vectorField mP(Md^fP);
    const vectorField mV(Md^fV);

    sumPatchMomentsP_ += sum(mP);
    sumPatchMomentsV_ += sum(mV);
    moment().boundaryFieldRef()[patchi] += mP + mV;
}


void Foam::functionObjects::forces::addToInternalField
(
    const labelList& cellIDs,
    const vectorField& Md,
    const vectorField& f
)
{
    auto& force = this->force();
    auto& moment = this->moment();

    forAll(cellIDs, i)
    {
        const label celli = cellIDs[i];

        sumInternalForces_ += f[i];
        force[celli] += f[i];

        const vector m(Md[i]^f[i]);
        sumInternalMoments_ += m;
        moment[celli] = m;
    }
}


void Foam::functionObjects::forces::createIntegratedDataFiles()
{
    if (!forceFilePtr_.valid())
    {
        forceFilePtr_ = createFile("force");
        writeIntegratedDataFileHeader("Force", forceFilePtr_());
    }

    if (!momentFilePtr_.valid())
    {
        momentFilePtr_ = createFile("moment");
        writeIntegratedDataFileHeader("Moment", momentFilePtr_());
    }
}


void Foam::functionObjects::forces::writeIntegratedDataFileHeader
(
    const word& header,
    OFstream& os
) const
{
    const auto& coordSys = coordSysPtr_();
    const auto vecDesc = [](const word& root)->string
    {
        return root + "_x " + root + "_y " + root + "_z";
    };
    writeHeader(os, header);
    writeHeaderValue(os, "CofR", coordSys.origin());
    writeHeader(os, "");
    writeCommented(os, "Time");
    writeTabbed(os, vecDesc("total"));
    writeTabbed(os, vecDesc("pressure"));
    writeTabbed(os, vecDesc("viscous"));

    if (porosity_)
    {
        writeTabbed(os, vecDesc("porous"));
    }

    os  << endl;
}


void Foam::functionObjects::forces::writeIntegratedDataFiles()
{
    const auto& coordSys = coordSysPtr_();

    writeIntegratedDataFile
    (
        coordSys.localVector(sumPatchForcesP_),
        coordSys.localVector(sumPatchForcesV_),
        coordSys.localVector(sumInternalForces_),
        forceFilePtr_()
    );

    writeIntegratedDataFile
    (
        coordSys.localVector(sumPatchMomentsP_),
        coordSys.localVector(sumPatchMomentsV_),
        coordSys.localVector(sumInternalMoments_),
        momentFilePtr_()
    );
}


void Foam::functionObjects::forces::writeIntegratedDataFile
(
    const vector& pres,
    const vector& vis,
    const vector& internal,
    OFstream& os
) const
{
    writeCurrentTime(os);

    writeValue(os, pres + vis + internal);
    writeValue(os, pres);
    writeValue(os, vis);

    if (porosity_)
    {
        writeValue(os, internal);
    }

    os  << endl;
}


void Foam::functionObjects::forces::logIntegratedData
(
    const string& descriptor,
    const vector& pres,
    const vector& vis,
    const vector& internal
) const
{
    if (!log)
    {
        return;
    }

    Log << "    Sum of " << descriptor.c_str() << nl
        << "        Total    : " << (pres + vis + internal) << nl
        << "        Pressure : " << pres << nl
        << "        Viscous  : " << vis << nl;

    if (porosity_)
    {
        Log << "        Porous   : " << internal << nl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::forces::forces
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    bool readFields
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name),
    sumPatchForcesP_(Zero),
    sumPatchForcesV_(Zero),
    sumPatchMomentsP_(Zero),
    sumPatchMomentsV_(Zero),
    sumInternalForces_(Zero),
    sumInternalMoments_(Zero),
    forceFilePtr_(),
    momentFilePtr_(),
    coordSysPtr_(nullptr),
    patchSet_(),
    rhoRef_(VGREAT),
    pRef_(0),
    pName_("p"),
    UName_("U"),
    rhoName_("rho"),
    fDName_("fD"),
    directForceDensity_(false),
    porosity_(false),
    writeFields_(false),
    initialised_(false)
{
    if (readFields)
    {
        read(dict);
        setCoordinateSystem(dict);
        Log << endl;
    }
}


Foam::functionObjects::forces::forces
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    bool readFields
)
:
    fvMeshFunctionObject(name, obr, dict),
    writeFile(mesh_, name),
    sumPatchForcesP_(Zero),
    sumPatchForcesV_(Zero),
    sumPatchMomentsP_(Zero),
    sumPatchMomentsV_(Zero),
    sumInternalForces_(Zero),
    sumInternalMoments_(Zero),
    forceFilePtr_(),
    momentFilePtr_(),
    coordSysPtr_(nullptr),
    patchSet_(),
    rhoRef_(VGREAT),
    pRef_(0),
    pName_("p"),
    UName_("U"),
    rhoName_("rho"),
    fDName_("fD"),
    directForceDensity_(false),
    porosity_(false),
    writeFields_(false),
    initialised_(false)
{
    if (readFields)
    {
        read(dict);
        setCoordinateSystem(dict);
        Log << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::forces::read(const dictionary& dict)
{
    if (!fvMeshFunctionObject::read(dict) || !writeFile::read(dict))
    {
        return false;
    }

    initialised_ = false;

    Info<< type() << " " << name() << ":" << endl;

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            dict.get<wordRes>("patches")
        );

    dict.readIfPresent("directForceDensity", directForceDensity_);
    if (directForceDensity_)
    {
        // Optional name entry for fD
        if (dict.readIfPresent<word>("fD", fDName_))
        {
            Info<< "    fD: " << fDName_ << endl;
        }
    }
    else
    {
        // Optional field name entries
        if (dict.readIfPresent<word>("p", pName_))
        {
            Info<< "    p: " << pName_ << endl;
        }
        if (dict.readIfPresent<word>("U", UName_))
        {
            Info<< "    U: " << UName_ << endl;
        }
        if (dict.readIfPresent<word>("rho", rhoName_))
        {
            Info<< "    rho: " << rhoName_ << endl;
        }

        // Reference density needed for incompressible calculations
        if (rhoName_ == "rhoInf")
        {
            rhoRef_ = dict.getCheck<scalar>("rhoInf", scalarMinMax::ge(SMALL));
            Info<< "    Freestream density (rhoInf) set to " << rhoRef_ << endl;
        }

        // Reference pressure, 0 by default
        if (dict.readIfPresent<scalar>("pRef", pRef_))
        {
            Info<< "    Reference pressure (pRef) set to " << pRef_ << endl;
        }
    }

    dict.readIfPresent("porosity", porosity_);
    if (porosity_)
    {
        Info<< "    Including porosity effects" << endl;
    }
    else
    {
        Info<< "    Not including porosity effects" << endl;
    }

    writeFields_ = dict.getOrDefault("writeFields", false);
    if (writeFields_)
    {
        Info<< "    Fields will be written" << endl;
    }


    return true;
}


void Foam::functionObjects::forces::calcForcesMoments()
{
    initialise();

    reset();

    const point& origin = coordSysPtr_->origin();

    if (directForceDensity_)
    {
        const auto& fD = lookupObject<volVectorField>(fDName_);

        const auto& Sfb = mesh_.Sf().boundaryField();

        for (const label patchi : patchSet_)
        {
            const vectorField& d = mesh_.C().boundaryField()[patchi];

            const vectorField Md(d - origin);

            const scalarField sA(mag(Sfb[patchi]));

            // Pressure force = surfaceUnitNormal*(surfaceNormal & forceDensity)
            const vectorField fP
            (
                Sfb[patchi]/sA
               *(
                    Sfb[patchi] & fD.boundaryField()[patchi]
                )
            );

            // Viscous force (total force minus pressure fP)
            const vectorField fV(sA*fD.boundaryField()[patchi] - fP);

            addToPatchFields(patchi, Md, fP, fV);
        }
    }
    else
    {
        const auto& p = lookupObject<volScalarField>(pName_);

        const auto& Sfb = mesh_.Sf().boundaryField();

        tmp<volSymmTensorField> tdevRhoReff = devRhoReff();
        const auto& devRhoReffb = tdevRhoReff().boundaryField();

        // Scale pRef by density for incompressible simulations
        const scalar rhoRef = rho(p);
        const scalar pRef = pRef_/rhoRef;

        for (const label patchi : patchSet_)
        {
            const vectorField& d = mesh_.C().boundaryField()[patchi];

            const vectorField Md(d - origin);

            const vectorField fP
            (
                rhoRef*Sfb[patchi]*(p.boundaryField()[patchi] - pRef)
            );

            const vectorField fV(Sfb[patchi] & devRhoReffb[patchi]);

            addToPatchFields(patchi, Md, fP, fV);
        }
    }

    if (porosity_)
    {
        const auto& U = lookupObject<volVectorField>(UName_);
        const volScalarField rho(this->rho());
        const volScalarField mu(this->mu());

        const auto models = obr_.lookupClass<porosityModel>();

        if (models.empty())
        {
            WarningInFunction
                << "Porosity effects requested, "
                << "but no porosity models found in the database"
                << endl;
        }

        forAllConstIters(models, iter)
        {
            // Non-const access required if mesh is changing
            auto& pm = const_cast<porosityModel&>(*iter());

            const vectorField fPTot(pm.force(U, rho, mu));

            const labelList& cellZoneIDs = pm.cellZoneIDs();

            for (const label zonei : cellZoneIDs)
            {
                const cellZone& cZone = mesh_.cellZones()[zonei];

                const vectorField d(mesh_.C(), cZone);
                const vectorField fP(fPTot, cZone);
                const vectorField Md(d - origin);

                addToInternalField(cZone, Md, fP);
            }
        }
    }

    reduce(sumPatchForcesP_, sumOp<vector>());
    reduce(sumPatchForcesV_, sumOp<vector>());
    reduce(sumPatchMomentsP_, sumOp<vector>());
    reduce(sumPatchMomentsV_, sumOp<vector>());
    reduce(sumInternalForces_, sumOp<vector>());
    reduce(sumInternalMoments_, sumOp<vector>());
}


Foam::vector Foam::functionObjects::forces::forceEff() const
{
    return sumPatchForcesP_ + sumPatchForcesV_ + sumInternalForces_;
}


Foam::vector Foam::functionObjects::forces::momentEff() const
{
    return sumPatchMomentsP_ + sumPatchMomentsV_ + sumInternalMoments_;
}


bool Foam::functionObjects::forces::execute()
{
    calcForcesMoments();

    Log << type() << " " << name() << " write:" << nl;

    const auto& coordSys = coordSysPtr_();

    const auto localFp(coordSys.localVector(sumPatchForcesP_));
    const auto localFv(coordSys.localVector(sumPatchForcesV_));
    const auto localFi(coordSys.localVector(sumInternalForces_));

    logIntegratedData("forces", localFp, localFv, localFi);

    const auto localMp(coordSys.localVector(sumPatchMomentsP_));
    const auto localMv(coordSys.localVector(sumPatchMomentsV_));
    const auto localMi(coordSys.localVector(sumInternalMoments_));

    logIntegratedData("moments", localMp, localMv, localMi);

    setResult("pressureForce", localFp);
    setResult("viscousForce", localFv);
    setResult("internalForce", localFi);
    setResult("pressureMoment", localMp);
    setResult("viscousMoment", localMv);
    setResult("internalMoment", localMi);

    return true;
}


bool Foam::functionObjects::forces::write()
{
    if (writeToFile())
    {
        Log << "    writing force and moment files." << endl;

        createIntegratedDataFiles();
        writeIntegratedDataFiles();
    }

    if (writeFields_)
    {
        Log << "    writing force and moment fields." << endl;

        force().write();
        moment().write();
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
