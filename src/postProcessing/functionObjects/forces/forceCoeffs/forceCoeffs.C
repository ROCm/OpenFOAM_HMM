/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd
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

#include "forceCoeffs.H"
#include "dictionary.H"
#include "Time.H"
#include "Pstream.H"
#include "IOmanip.H"
#include "fvMesh.H"
#include "dimensionedTypes.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(forceCoeffs, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::forceCoeffs::createFiles()
{
    // Note: Only possible to create bin files after bins have been initialised

    if (writeToFile() && !coeffFilePtr_.valid())
    {
        coeffFilePtr_ = createFile("coefficient");
        writeIntegratedHeader("Coefficients", coeffFilePtr_());

        if (nBin_ > 1)
        {
            CmBinFilePtr_ = createFile("CmBin");
            writeBinHeader("Moment coefficient bins", CmBinFilePtr_());
            CdBinFilePtr_ = createFile("CdBin");
            writeBinHeader("Drag coefficient bins", CdBinFilePtr_());
            ClBinFilePtr_ = createFile("ClBin");
            writeBinHeader("Lift coefficient bins", ClBinFilePtr_());
        }
    }
}


void Foam::forceCoeffs::writeIntegratedHeader
(
    const word& header,
    Ostream& os
) const
{
    writeHeader(os, "Force coefficients");
    writeHeaderValue(os, "liftDir", liftDir_);
    writeHeaderValue(os, "dragDir", dragDir_);
    writeHeaderValue(os, "pitchAxis", pitchAxis_);
    writeHeaderValue(os, "magUInf", magUInf_);
    writeHeaderValue(os, "lRef", lRef_);
    writeHeaderValue(os, "Aref", Aref_);
    writeHeaderValue(os, "CofR", coordSys_.origin());
    writeHeader(os, "");
    writeCommented(os, "Time");
    writeTabbed(os, "Cm");
    writeTabbed(os, "Cd");
    writeTabbed(os, "Cl");
    writeTabbed(os, "Cl(f)");
    writeTabbed(os, "Cl(r)");
    os  << endl;
}


void Foam::forceCoeffs::writeBinHeader
(
    const word& header,
    Ostream& os
) const
{
    writeHeader(os, header);
    writeHeaderValue(os, "bins", nBin_);
    writeHeaderValue(os, "start", binMin_);
    writeHeaderValue(os, "delta", binDx_);
    writeHeaderValue(os, "direction", binDir_);

    vectorField binPoints(nBin_);
    writeCommented(os, "x co-ords  :");
    forAll(binPoints, pointI)
    {
        binPoints[pointI] = (binMin_ + (pointI + 1)*binDx_)*binDir_;
        os << tab << binPoints[pointI].x();
    }
    os << nl;

    writeCommented(os, "y co-ords  :");
    forAll(binPoints, pointI)
    {
        os << tab << binPoints[pointI].y();
    }
    os << nl;

    writeCommented(os, "z co-ords  :");
    forAll(binPoints, pointI)
    {
        os << tab << binPoints[pointI].z();
    }
    os << nl;

    writeHeader(os, "");
    writeCommented(os, "Time");

    for (label j = 0; j < nBin_; j++)
    {
        word jn(Foam::name(j) + ':');
        writeTabbed(os, jn + "total");
        writeTabbed(os, jn + "pressure");
        writeTabbed(os, jn + "viscous");

        if (porosity_)
        {
            writeTabbed(os, jn + "porous");
        }
    }

    os  << endl;
}


void Foam::forceCoeffs::writeIntegratedData
(
    const word& title,
    const List<Field<scalar> >& coeff
) const
{
    scalar pressure = sum(coeff[0]);
    scalar viscous = sum(coeff[1]);
    scalar porous = sum(coeff[2]);
    scalar total = pressure + viscous + porous;

    if (log_)
    {
        Info<< "        " << title << "       : " << total << token::TAB
            << "("
            << "pressure: " << pressure << token::TAB
            << "viscous: " << viscous;

        if (porosity_)
        {
            Info<< token::TAB << "porous: " << porous;
        }

        Info<< ")" << endl;
    }
}


void Foam::forceCoeffs::writeBinData
(
    const List<Field<scalar> > coeffs,
    Ostream& os
) const
{
    os  << obr_.time().value();

    for (label binI = 0; binI < nBin_; binI++)
    {
        scalar total = coeffs[0][binI] + coeffs[1][binI] + coeffs[2][binI];

        os  << tab << total << tab << coeffs[0][binI] << tab << coeffs[1][binI];

        if (porosity_)
        {
            os  << tab << coeffs[2][binI];
        }
    }

    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::forceCoeffs::forceCoeffs
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles,
    const bool readFields
)
:
    forces(name, obr, dict, loadFromFiles, false),
    liftDir_(vector::zero),
    dragDir_(vector::zero),
    pitchAxis_(vector::zero),
    magUInf_(0.0),
    lRef_(0.0),
    Aref_(0.0),
    coeffFilePtr_(),
    CmBinFilePtr_(),
    CdBinFilePtr_(),
    ClBinFilePtr_()
{
    if (readFields)
    {
        read(dict);
        if (log_) Info << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::forceCoeffs::~forceCoeffs()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::forceCoeffs::read(const dictionary& dict)
{
    if (!active_)
    {
        return;
    }

    forces::read(dict);

    // Directions for lift and drag forces, and pitch moment
    dict.lookup("liftDir") >> liftDir_;
    dict.lookup("dragDir") >> dragDir_;
    dict.lookup("pitchAxis") >> pitchAxis_;

    // Free stream velocity magnitude
    dict.lookup("magUInf") >> magUInf_;

    // Reference length and area scales
    dict.lookup("lRef") >> lRef_;
    dict.lookup("Aref") >> Aref_;

    if (writeFields_)
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        tmp<volVectorField> tforceCoeff
        (
            new volVectorField
            (
                IOobject
                (
                    fieldName("forceCoeff"),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedVector("0", dimless, vector::zero)
            )
        );

        obr_.store(tforceCoeff.ptr());

        tmp<volVectorField> tmomentCoeff
        (
            new volVectorField
            (
                IOobject
                (
                    fieldName("momentCoeff"),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedVector("0", dimless, vector::zero)
            )
        );

        obr_.store(tmomentCoeff.ptr());
    }
}


void Foam::forceCoeffs::execute()
{
    if (!active_)
    {
        return;
    }

    forces::calcForcesMoment();

    createFiles();

    scalar pDyn = 0.5*rhoRef_*magUInf_*magUInf_;

    // Storage for pressure, viscous and porous contributions to coeffs
    List<Field<scalar> > momentCoeffs(3);
    List<Field<scalar> > dragCoeffs(3);
    List<Field<scalar> > liftCoeffs(3);
    forAll(liftCoeffs, i)
    {
        momentCoeffs[i].setSize(nBin_);
        dragCoeffs[i].setSize(nBin_);
        liftCoeffs[i].setSize(nBin_);
    }

    // Calculate coefficients
    scalar CmTot = 0;
    scalar CdTot = 0;
    scalar ClTot = 0;
    forAll(liftCoeffs, i)
    {
        momentCoeffs[i] = (moment_[i] & pitchAxis_)/(Aref_*pDyn*lRef_);
        dragCoeffs[i] = (force_[i] & dragDir_)/(Aref_*pDyn);
        liftCoeffs[i] = (force_[i] & liftDir_)/(Aref_*pDyn);

        CmTot += sum(momentCoeffs[i]);
        CdTot += sum(dragCoeffs[i]);
        ClTot += sum(liftCoeffs[i]);
    }

    scalar ClfTot = ClTot/2.0 + CmTot;
    scalar ClrTot = ClTot/2.0 - CmTot;

    if (log_) Info
        << type() << " " << name_ << " output:" << nl
        << "    Coefficients" << nl;

    writeIntegratedData("Cm", momentCoeffs);
    writeIntegratedData("Cd", dragCoeffs);
    writeIntegratedData("Cl", liftCoeffs);

    if (log_) Info
        << "        Cl(f)    : " << ClfTot << nl
        << "        Cl(r)    : " << ClrTot << nl
        << endl;

    if (writeToFile())
    {
        writeTime(coeffFilePtr_());
        coeffFilePtr_()
            << tab << CmTot << tab  << CdTot
            << tab << ClTot << tab << ClfTot << tab << ClrTot << endl;


        if (nBin_ > 1)
        {
            if (binCumulative_)
            {
                forAll(liftCoeffs, i)
                {
                    for (label binI = 1; binI < nBin_; binI++)
                    {
                        liftCoeffs[i][binI] += liftCoeffs[i][binI-1];
                        dragCoeffs[i][binI] += dragCoeffs[i][binI-1];
                        momentCoeffs[i][binI] += momentCoeffs[i][binI-1];
                    }
                }
            }

            writeBinData(dragCoeffs, CdBinFilePtr_());
            writeBinData(liftCoeffs, ClBinFilePtr_());
            writeBinData(momentCoeffs, CmBinFilePtr_());
        }
    }

    // Write state/results information
    {
        setResult("Cm", CmTot);
        setResult("Cd", CdTot);
        setResult("Cl", ClTot);
        setResult("Cl(f)", ClfTot);
        setResult("Cl(r)", ClrTot);
    }

    if (writeFields_)
    {
        const volVectorField& force =
            obr_.lookupObject<volVectorField>(fieldName("force"));

        const volVectorField& moment =
            obr_.lookupObject<volVectorField>(fieldName("moment"));

        volVectorField& forceCoeff =
            const_cast<volVectorField&>
            (
                obr_.lookupObject<volVectorField>(fieldName("forceCoeff"))
            );

        volVectorField& momentCoeff =
            const_cast<volVectorField&>
            (
                obr_.lookupObject<volVectorField>(fieldName("momentCoeff"))
            );

        dimensionedScalar f0("f0", dimForce, Aref_*pDyn);
        dimensionedScalar m0("m0", dimForce*dimLength, Aref_*lRef_*pDyn);

        forceCoeff == force/f0;
        momentCoeff == moment/m0;
    }
}


void Foam::forceCoeffs::end()
{
    if (active_)
    {
        execute();
    }
}


void Foam::forceCoeffs::timeSet()
{
    // Do nothing
}


void Foam::forceCoeffs::write()
{
    if (!active_)
    {
        return;
    }

    if (writeFields_)
    {
        const volVectorField& forceCoeff =
            obr_.lookupObject<volVectorField>(fieldName("forceCoeff"));

        const volVectorField& momentCoeff =
            obr_.lookupObject<volVectorField>(fieldName("momentCoeff"));

        forceCoeff.write();
        momentCoeff.write();
    }
}


// ************************************************************************* //
