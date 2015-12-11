/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "wordReList.H"
#include "fvcGrad.H"
#include "porosityModel.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(forces, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::word Foam::forces::fieldName(const word& name) const
{
    return name_ + ":" + name;
}


void Foam::forces::createFiles()
{
    // Note: Only possible to create bin files after bins have been initialised

    if (writeToFile() && !forceFilePtr_.valid())
    {
        forceFilePtr_ = createFile("force");
        writeIntegratedHeader("Force", forceFilePtr_());
        momentFilePtr_ = createFile("moment");
        writeIntegratedHeader("Moment", momentFilePtr_());

        if (nBin_ > 1)
        {
            forceBinFilePtr_ = createFile("forceBin");
            writeBinHeader("Force", forceBinFilePtr_());
            momentBinFilePtr_ = createFile("momentBin");
            writeBinHeader("Moment", momentBinFilePtr_());
        }

        if (localSystem_)
        {
            localForceFilePtr_ = createFile("localForce");
            writeIntegratedHeader("Force", localForceFilePtr_());
            localMomentFilePtr_ = createFile("localMoment");
            writeIntegratedHeader("Moment", localMomentFilePtr_());

            if (nBin_ > 1)
            {
                localForceBinFilePtr_ = createFile("localForceBin");
                writeBinHeader("Force", localForceBinFilePtr_());
                localMomentBinFilePtr_ = createFile("localMomentBin");
                writeBinHeader("Moment", localMomentBinFilePtr_());
            }
        }
    }
}


void Foam::forces::writeIntegratedHeader
(
    const word& header,
    Ostream& os
) const
{
    writeHeader(os, header);
    writeHeaderValue(os, "CofR", coordSys_.origin());
    writeHeader(os, "");
    writeCommented(os, "Time");
    writeTabbed(os, "(total_x total_y total_z)");
    writeTabbed(os, "(pressure_x pressure_y pressure_z)");
    writeTabbed(os, "(viscous_x viscous_y viscous_z)");

    if (porosity_)
    {
        writeTabbed(os, "(porous_x porous_y porous_z)");
    }

    os  << endl;
}


void Foam::forces::writeBinHeader(const word& header, Ostream& os) const
{
    writeHeader(os, header + " bins");
    writeHeaderValue(os, "bins", nBin_);
    writeHeaderValue(os, "start", binMin_);
    writeHeaderValue(os, "delta", binDx_);
    writeHeaderValue(os, "direction", binDir_);

    vectorField binPoints(nBin_);
    writeCommented(os, "x co-ords  :");
    forAll(binPoints, pointI)
    {
        binPoints[pointI] = (binMin_ + (pointI + 1)*binDx_)*binDir_;
        os  << tab << binPoints[pointI].x();
    }
    os  << nl;

    writeCommented(os, "y co-ords  :");
    forAll(binPoints, pointI)
    {
        os  << tab << binPoints[pointI].y();
    }
    os  << nl;

    writeCommented(os, "z co-ords  :");
    forAll(binPoints, pointI)
    {
        os  << tab << binPoints[pointI].z();
    }
    os  << nl;

    writeHeader(os, "");
    writeCommented(os, "Time");

    for (label j = 0; j < nBin_; j++)
    {
        const word jn(Foam::name(j) + ':');
        os  << tab << jn << "(total_x total_y total_z)"
            << tab << jn << "(pressure_x pressure_y pressure_z)"
            << tab << jn << "(viscous_x viscous_y viscous_z)";

        if (porosity_)
        {
            os  << tab << jn << "(porous_x porous_y porous_z)";
        }
    }

    os << endl;
}


void Foam::forces::initialise()
{
    if (initialised_ || !active_)
    {
        return;
    }

    if (directForceDensity_)
    {
        if (!obr_.foundObject<volVectorField>(fDName_))
        {
            active_ = false;
            WarningInFunction
                << "Could not find " << fDName_ << " in database." << nl
                << "    De-activating forces."
                << endl;
        }
    }
    else
    {
        if
        (
            !obr_.foundObject<volVectorField>(UName_)
         || !obr_.foundObject<volScalarField>(pName_)
         || (
                rhoName_ != "rhoInf"
             && !obr_.foundObject<volScalarField>(rhoName_)
            )
        )
        {
            active_ = false;

            WarningInFunction
                << "Could not find " << UName_ << ", " << pName_;

            if (rhoName_ != "rhoInf")
            {
                Info<< " or " << rhoName_;
            }

            Info<< " in database." << nl
                << "    De-activating forces." << endl;
        }
    }

    initialiseBins();

    initialised_ = true;
}


void Foam::forces::initialiseBins()
{
    if (!active_)
    {
        return;
    }

    if (nBin_ > 1)
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();

        // Determine extents of patches
        binMin_ = GREAT;
        scalar binMax = -GREAT;
        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchI = iter.key();
            const polyPatch& pp = pbm[patchI];
            scalarField d(pp.faceCentres() & binDir_);
            binMin_ = min(min(d), binMin_);
            binMax = max(max(d), binMax);
        }

        // Include porosity
        if (porosity_)
        {
            const HashTable<const porosityModel*> models =
                obr_.lookupClass<porosityModel>();

            const scalarField dd(mesh.C() & binDir_);

            forAllConstIter(HashTable<const porosityModel*>, models, iter)
            {
                const porosityModel& pm = *iter();
                const labelList& cellZoneIDs = pm.cellZoneIDs();

                forAll(cellZoneIDs, i)
                {
                    label zoneI = cellZoneIDs[i];
                    const cellZone& cZone = mesh.cellZones()[zoneI];
                    const scalarField d(dd, cZone);
                    binMin_ = min(min(d), binMin_);
                    binMax = max(max(d), binMax);
                }
            }
        }

        reduce(binMin_, minOp<scalar>());
        reduce(binMax, maxOp<scalar>());

        // Slightly boost binMax so that region of interest is fully
        // within bounds
        binMax = 1.0001*(binMax - binMin_) + binMin_;

        binDx_ = (binMax - binMin_)/scalar(nBin_);

        // Create the bin points used for writing
        binPoints_.setSize(nBin_);
        forAll(binPoints_, i)
        {
            binPoints_[i] = (i + 0.5)*binDir_*binDx_;
        }

        // Allocate storage for forces and moments
        forAll(force_, i)
        {
            force_[i].setSize(nBin_);
            moment_[i].setSize(nBin_);
        }
    }
}


void Foam::forces::resetFields()
{
    force_[0] = vector::zero;
    force_[1] = vector::zero;
    force_[2] = vector::zero;

    moment_[0] = vector::zero;
    moment_[1] = vector::zero;
    moment_[2] = vector::zero;

    if (writeFields_)
    {
        volVectorField& force =
            const_cast<volVectorField&>
            (
                obr_.lookupObject<volVectorField>(fieldName("force"))
            );

        force == dimensionedVector("0", force.dimensions(), vector::zero);

        volVectorField& moment =
            const_cast<volVectorField&>
            (
                obr_.lookupObject<volVectorField>(fieldName("moment"))
            );

        moment == dimensionedVector("0", moment.dimensions(), vector::zero);
    }
}


Foam::tmp<Foam::volSymmTensorField> Foam::forces::devRhoReff() const
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    if (obr_.foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        const cmpTurbModel& turb =
            obr_.lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);

        return turb.devRhoReff();
    }
    else if (obr_.foundObject<icoTurbModel>(icoTurbModel::propertiesName))
    {
        const incompressible::turbulenceModel& turb =
            obr_.lookupObject<icoTurbModel>(icoTurbModel::propertiesName);

        return rho()*turb.devReff();
    }
    else if (obr_.foundObject<fluidThermo>(fluidThermo::typeName))
    {
        const fluidThermo& thermo =
            obr_.lookupObject<fluidThermo>(fluidThermo::typeName);

        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

        return -thermo.mu()*dev(twoSymm(fvc::grad(U)));
    }
    else if
    (
        obr_.foundObject<transportModel>("transportProperties")
    )
    {
        const transportModel& laminarT =
            obr_.lookupObject<transportModel>("transportProperties");

        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

        return -rho()*laminarT.nu()*dev(twoSymm(fvc::grad(U)));
    }
    else if (obr_.foundObject<dictionary>("transportProperties"))
    {
        const dictionary& transportProperties =
             obr_.lookupObject<dictionary>("transportProperties");

        dimensionedScalar nu(transportProperties.lookup("nu"));

        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

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


Foam::tmp<Foam::volScalarField> Foam::forces::mu() const
{
    if (obr_.foundObject<fluidThermo>(basicThermo::dictName))
    {
        const fluidThermo& thermo =
             obr_.lookupObject<fluidThermo>(basicThermo::dictName);

        return thermo.mu();
    }
    else if
    (
        obr_.foundObject<transportModel>("transportProperties")
    )
    {
        const transportModel& laminarT =
            obr_.lookupObject<transportModel>("transportProperties");

        return rho()*laminarT.nu();
    }
    else if (obr_.foundObject<dictionary>("transportProperties"))
    {
        const dictionary& transportProperties =
             obr_.lookupObject<dictionary>("transportProperties");

        dimensionedScalar nu(transportProperties.lookup("nu"));

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


Foam::tmp<Foam::volScalarField> Foam::forces::rho() const
{
    if (rhoName_ == "rhoInf")
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar("rho", dimDensity, rhoRef_)
            )
        );
    }
    else
    {
        return(obr_.lookupObject<volScalarField>(rhoName_));
    }
}


Foam::scalar Foam::forces::rho(const volScalarField& p) const
{
    if (p.dimensions() == dimPressure)
    {
        return 1.0;
    }
    else
    {
        if (rhoName_ != "rhoInf")
        {
            FatalErrorInFunction
                << "Dynamic pressure is expected but kinematic is provided."
                << exit(FatalError);
        }

        return rhoRef_;
    }
}


void Foam::forces::applyBins
(
    const vectorField& Md,
    const vectorField& fN,
    const vectorField& fT,
    const vectorField& fP,
    const vectorField& d
)
{
    if (nBin_ == 1)
    {
        force_[0][0] += sum(fN);
        force_[1][0] += sum(fT);
        force_[2][0] += sum(fP);
        moment_[0][0] += sum(Md^fN);
        moment_[1][0] += sum(Md^fT);
        moment_[2][0] += sum(Md^fP);
    }
    else
    {
        scalarField dd((d & binDir_) - binMin_);

        forAll(dd, i)
        {
            label bini = min(max(floor(dd[i]/binDx_), 0), force_[0].size() - 1);

            force_[0][bini] += fN[i];
            force_[1][bini] += fT[i];
            force_[2][bini] += fP[i];
            moment_[0][bini] += Md[i]^fN[i];
            moment_[1][bini] += Md[i]^fT[i];
            moment_[2][bini] += Md[i]^fP[i];
        }
    }
}


void Foam::forces::addToFields
(
    const label patchI,
    const vectorField& Md,
    const vectorField& fN,
    const vectorField& fT,
    const vectorField& fP
)
{
    if (!writeFields_)
    {
        return;
    }

    volVectorField& force =
        const_cast<volVectorField&>
        (
            obr_.lookupObject<volVectorField>(fieldName("force"))
        );

    vectorField& pf = force.boundaryField()[patchI];
    pf += fN + fT + fP;

    volVectorField& moment =
        const_cast<volVectorField&>
        (
            obr_.lookupObject<volVectorField>(fieldName("moment"))
        );

    vectorField& pm = moment.boundaryField()[patchI];
    pm += Md;
}


void Foam::forces::addToFields
(
    const labelList& cellIDs,
    const vectorField& Md,
    const vectorField& fN,
    const vectorField& fT,
    const vectorField& fP
)
{
    if (!writeFields_)
    {
        return;
    }

    volVectorField& force =
        const_cast<volVectorField&>
        (
            obr_.lookupObject<volVectorField>(fieldName("force"))
        );

    volVectorField& moment =
        const_cast<volVectorField&>
        (
            obr_.lookupObject<volVectorField>(fieldName("moment"))
        );

    forAll(cellIDs, i)
    {
        label cellI = cellIDs[i];
        force[cellI] += fN[i] + fT[i] + fP[i];
        moment[cellI] += Md[i];
    }
}


void Foam::forces::writeIntegratedForceMoment
(
    const string& descriptor,
    const vectorField& fm0,
    const vectorField& fm1,
    const vectorField& fm2,
    autoPtr<OFstream>& osPtr
) const
{
    vector pressure = sum(fm0);
    vector viscous = sum(fm1);
    vector porous = sum(fm2);
    vector total = pressure + viscous + porous;

    if (log_)
    {
        Info<< "    Sum of " << descriptor.c_str() << nl
            << "        Total    : " << total << nl
            << "        Pressure : " << pressure << nl
            << "        Viscous  : " << viscous << nl;

        if (porosity_)
        {
            Info<< "        Porous   : " << porous << nl;
        }
    }

    if (writeToFile())
    {
        Ostream& os = osPtr();

        os  << obr_.time().value()
            << tab << total
            << tab << pressure
            << tab << viscous;

        if (porosity_)
        {
            os  << tab << porous;
        }

        os  << endl;
    }
}


void Foam::forces::writeForces()
{
    if (log_) Info << type() << " " << name_ << " output:" << nl;

    writeIntegratedForceMoment
    (
        "forces",
        force_[0],
        force_[1],
        force_[2],
        forceFilePtr_
    );

    writeIntegratedForceMoment
    (
        "moments",
        moment_[0],
        moment_[1],
        moment_[2],
        momentFilePtr_
    );

    if (localSystem_)
    {
        writeIntegratedForceMoment
        (
            "local forces",
            coordSys_.localVector(force_[0]),
            coordSys_.localVector(force_[1]),
            coordSys_.localVector(force_[2]),
            localForceFilePtr_
        );

        writeIntegratedForceMoment
        (
            "local moments",
            coordSys_.localVector(moment_[0]),
            coordSys_.localVector(moment_[1]),
            coordSys_.localVector(moment_[2]),
            localMomentFilePtr_
        );
    }

    if (log_) Info << endl;
}


void Foam::forces::writeBinnedForceMoment
(
    const List<Field<vector> >& fm,
    autoPtr<OFstream>& osPtr
) const
{
    if ((nBin_ == 1) || !writeToFile())
    {
        return;
    }

    List<Field<vector> > f(fm);

    if (binCumulative_)
    {
        for (label i = 1; i < f[0].size(); i++)
        {
            f[0][i] += f[0][i-1];
            f[1][i] += f[1][i-1];
            f[2][i] += f[2][i-1];
        }
    }

    Ostream& os = osPtr();

    writeTime(os);

    forAll(f[0], i)
    {
        vector total = f[0][i] + f[1][i] + f[2][i];

        os  << tab << total
            << tab << f[0][i]
            << tab << f[1][i];

        if (porosity_)
        {
            os  << tab << f[2][i];
        }
    }

    os  << nl;
}


void Foam::forces::writeBins()
{
    writeBinnedForceMoment(force_, forceBinFilePtr_);
    writeBinnedForceMoment(moment_, momentBinFilePtr_);

    if (localSystem_)
    {
        List<Field<vector> > lf(3);
        List<Field<vector> > lm(3);
        lf[0] = coordSys_.localVector(force_[0]);
        lf[1] = coordSys_.localVector(force_[1]);
        lf[2] = coordSys_.localVector(force_[2]);
        lm[0] = coordSys_.localVector(moment_[0]);
        lm[1] = coordSys_.localVector(moment_[1]);
        lm[2] = coordSys_.localVector(moment_[2]);

        writeBinnedForceMoment(lf, localForceBinFilePtr_);
        writeBinnedForceMoment(lm, localMomentBinFilePtr_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::forces::forces
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles,
    const bool readFields
)
:
    functionObjectState(obr, name),
    functionObjectFile(obr, name),
    obr_(obr),
    log_(true),
    force_(3),
    moment_(3),
    forceFilePtr_(),
    momentFilePtr_(),
    forceBinFilePtr_(),
    momentBinFilePtr_(),
    localForceFilePtr_(),
    localMomentFilePtr_(),
    localForceBinFilePtr_(),
    localMomentBinFilePtr_(),
    patchSet_(),
    pName_(word::null),
    UName_(word::null),
    rhoName_(word::null),
    directForceDensity_(false),
    fDName_(""),
    rhoRef_(VGREAT),
    pRef_(0),
    coordSys_(),
    localSystem_(false),
    porosity_(false),
    nBin_(1),
    binDir_(vector::zero),
    binDx_(0.0),
    binMin_(GREAT),
    binPoints_(),
    binCumulative_(true),
    writeFields_(false),
    initialised_(false)
{
    // Check if the available mesh is an fvMesh otherise deactivate
    if (setActive<fvMesh>())
    {
        if (readFields)
        {
            read(dict);
            if (log_) Info << endl;
        }
    }
}


Foam::forces::forces
(
    const word& name,
    const objectRegistry& obr,
    const labelHashSet& patchSet,
    const word& pName,
    const word& UName,
    const word& rhoName,
    const scalar rhoInf,
    const scalar pRef,
    const coordinateSystem& coordSys
)
:
    functionObjectState(obr, name),
    functionObjectFile(obr, name),
    obr_(obr),
    log_(true),
    force_(3),
    moment_(3),
    forceFilePtr_(),
    momentFilePtr_(),
    forceBinFilePtr_(),
    momentBinFilePtr_(),
    localForceFilePtr_(),
    localMomentFilePtr_(),
    localForceBinFilePtr_(),
    localMomentBinFilePtr_(),
    patchSet_(patchSet),
    pName_(pName),
    UName_(UName),
    rhoName_(rhoName),
    directForceDensity_(false),
    fDName_(""),
    rhoRef_(rhoInf),
    pRef_(pRef),
    coordSys_(coordSys),
    localSystem_(false),
    porosity_(false),
    nBin_(1),
    binDir_(vector::zero),
    binDx_(0.0),
    binMin_(GREAT),
    binPoints_(),
    binCumulative_(true),
    writeFields_(false),
    initialised_(false)
{
    // Turn off writing to file
    writeToFile_ = false;

    forAll(force_, i)
    {
        force_[i].setSize(nBin_);
        moment_[i].setSize(nBin_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::forces::~forces()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::forces::read(const dictionary& dict)
{
    if (!active_)
    {
        return;
    }

    functionObjectFile::read(dict);

    initialised_ = false;

    log_ = dict.lookupOrDefault<Switch>("log", false);

    if (log_) Info << type() << " " << name_ << ":" << nl;

    directForceDensity_ = dict.lookupOrDefault("directForceDensity", false);

    const fvMesh& mesh = refCast<const fvMesh>(obr_);
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    patchSet_ = pbm.patchSet(wordReList(dict.lookup("patches")));

    if (directForceDensity_)
    {
        // Optional entry for fDName
        fDName_ = dict.lookupOrDefault<word>("fDName", "fD");
    }
    else
    {
        // Optional entries U and p
        pName_ = dict.lookupOrDefault<word>("pName", "p");
        UName_ = dict.lookupOrDefault<word>("UName", "U");
        rhoName_ = dict.lookupOrDefault<word>("rhoName", "rho");

        // Reference density needed for incompressible calculations
        rhoRef_ = readScalar(dict.lookup("rhoInf"));

        // Reference pressure, 0 by default
        pRef_ = dict.lookupOrDefault<scalar>("pRef", 0.0);
    }

    coordSys_.clear();

    // Centre of rotation for moment calculations
    // specified directly, from coordinate system, or implicitly (0 0 0)
    if (!dict.readIfPresent<point>("CofR", coordSys_.origin()))
    {
        coordSys_ = coordinateSystem(obr_, dict);
        localSystem_ = true;
    }

    dict.readIfPresent("porosity", porosity_);
    if (porosity_)
    {
        if (log_) Info << "    Including porosity effects" << endl;
    }
    else
    {
        if (log_) Info << "    Not including porosity effects" << endl;
    }

    if (dict.found("binData"))
    {
        const dictionary& binDict(dict.subDict("binData"));
        binDict.lookup("nBin") >> nBin_;

        if (nBin_ < 0)
        {
            FatalIOErrorInFunction(dict)
                << "Number of bins (nBin) must be zero or greater"
                << exit(FatalIOError);
        }
        else if (nBin_ == 0)
        {
            nBin_ = 1;
        }
        else
        {
            binDict.lookup("cumulative") >> binCumulative_;
            binDict.lookup("direction") >> binDir_;
            binDir_ /= mag(binDir_);
        }
    }

    if (nBin_ == 1)
    {
        // Allocate storage for forces and moments
        forAll(force_, i)
        {
            force_[i].setSize(1);
            moment_[i].setSize(1);
        }
    }

    writeFields_ = dict.lookupOrDefault("writeFields", false);

    if (writeFields_)
    {
        if (log_) Info << "    Fields will be written" << endl;

        tmp<volVectorField> tforce
        (
            new volVectorField
            (
                IOobject
                (
                    fieldName("force"),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedVector("0", dimForce, vector::zero)
            )
        );

        obr_.store(tforce.ptr());

        tmp<volVectorField> tmoment
        (
            new volVectorField
            (
                IOobject
                (
                    fieldName("moment"),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedVector("0", dimForce*dimLength, vector::zero)
            )
        );

        obr_.store(tmoment.ptr());
    }
}


void Foam::forces::execute()
{
    if (!active_)
    {
        return;
    }

    // calcForcesMoment may have reset the active flag - need to re-check
    calcForcesMoment();

    if (!active_)
    {
        return;
    }

    if (Pstream::master())
    {
        createFiles();

        writeForces();

        writeBins();

        if (log_) Info << endl;
    }

    // write state/results information
    setResult("normalForce", sum(force_[0]));
    setResult("tangentialForce", sum(force_[1]));
    setResult("porousForce", sum(force_[2]));

    setResult("normalMoment", sum(moment_[0]));
    setResult("tangentialMoment", sum(moment_[1]));
    setResult("porousMoment", sum(moment_[2]));
}


void Foam::forces::end()
{
    // Do nothing
}


void Foam::forces::timeSet()
{
    // Do nothing
}


void Foam::forces::write()
{
    if (!active_)
    {
        return;
    }

    if (writeFields_)
    {
        obr_.lookupObject<volVectorField>(fieldName("force")).write();
        obr_.lookupObject<volVectorField>(fieldName("moment")).write();
    }
}


void Foam::forces::calcForcesMoment()
{
    if (!active_)
    {
        return;
    }

    initialise();

    // Initialise may have reset the active flag - need to re-check
    if (!active_)
    {
        return;
    }

    resetFields();

    if (directForceDensity_)
    {
        const volVectorField& fD = obr_.lookupObject<volVectorField>(fDName_);

        const fvMesh& mesh = fD.mesh();

        const surfaceVectorField::GeometricBoundaryField& Sfb =
            mesh.Sf().boundaryField();

        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchI = iter.key();

            vectorField Md
            (
                mesh.C().boundaryField()[patchI] - coordSys_.origin()
            );

            scalarField sA(mag(Sfb[patchI]));

            // Normal force = surfaceUnitNormal*(surfaceNormal & forceDensity)
            vectorField fN
            (
                Sfb[patchI]/sA
               *(
                    Sfb[patchI] & fD.boundaryField()[patchI]
                )
            );

            // Tangential force (total force minus normal fN)
            vectorField fT(sA*fD.boundaryField()[patchI] - fN);

            //- Porous force
            vectorField fP(Md.size(), vector::zero);

            addToFields(patchI, Md, fN, fT, fP);

            applyBins(Md, fN, fT, fP, mesh.C().boundaryField()[patchI]);
        }
    }
    else
    {
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
        const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);

        const fvMesh& mesh = U.mesh();

        const surfaceVectorField::GeometricBoundaryField& Sfb =
            mesh.Sf().boundaryField();

        tmp<volSymmTensorField> tdevRhoReff = devRhoReff();
        const volSymmTensorField::GeometricBoundaryField& devRhoReffb
            = tdevRhoReff().boundaryField();

        // Scale pRef by density for incompressible simulations
        scalar pRef = pRef_/rho(p);

        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchI = iter.key();

            vectorField Md
            (
                mesh.C().boundaryField()[patchI] - coordSys_.origin()
            );

            vectorField fN
            (
                rho(p)*Sfb[patchI]*(p.boundaryField()[patchI] - pRef)
            );

            vectorField fT(Sfb[patchI] & devRhoReffb[patchI]);

            vectorField fP(Md.size(), vector::zero);

            addToFields(patchI, Md, fN, fT, fP);

            applyBins(Md, fN, fT, fP, mesh.C().boundaryField()[patchI]);
        }
    }

    if (porosity_)
    {
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
        const volScalarField rho(this->rho());
        const volScalarField mu(this->mu());

        const fvMesh& mesh = U.mesh();

        const HashTable<const porosityModel*> models =
            obr_.lookupClass<porosityModel>();

        if (models.empty())
        {
            WarningInFunction
                << "Porosity effects requested, but no porosity models found "
                << "in the database"
                << endl;
        }

        forAllConstIter(HashTable<const porosityModel*>, models, iter)
        {
            // non-const access required if mesh is changing
            porosityModel& pm = const_cast<porosityModel&>(*iter());

            vectorField fPTot(pm.force(U, rho, mu));

            const labelList& cellZoneIDs = pm.cellZoneIDs();

            forAll(cellZoneIDs, i)
            {
                label zoneI = cellZoneIDs[i];
                const cellZone& cZone = mesh.cellZones()[zoneI];

                const vectorField d(mesh.C(), cZone);
                const vectorField fP(fPTot, cZone);
                const vectorField Md(d - coordSys_.origin());

                const vectorField fDummy(Md.size(), vector::zero);

                addToFields(cZone, Md, fDummy, fDummy, fP);

                applyBins(Md, fDummy, fDummy, fP, d);
            }
        }
    }

    Pstream::listCombineGather(force_, plusEqOp<vectorField>());
    Pstream::listCombineGather(moment_, plusEqOp<vectorField>());
    Pstream::listCombineScatter(force_);
    Pstream::listCombineScatter(moment_);
}


Foam::vector Foam::forces::forceEff() const
{
    return sum(force_[0]) + sum(force_[1]) + sum(force_[2]);
}


Foam::vector Foam::forces::momentEff() const
{
    return sum(moment_[0]) + sum(moment_[1]) + sum(moment_[2]);
}


// ************************************************************************* //
