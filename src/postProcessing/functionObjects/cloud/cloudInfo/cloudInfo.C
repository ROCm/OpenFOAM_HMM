/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2015 OpenFOAM Foundation
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

#include "cloudInfo.H"
#include "dictionary.H"
#include "kinematicCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cloudInfo, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::cloudInfo::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "Cloud information");
    writeCommented(os, "Time");
    writeTabbed(os, "nParcels");
    writeTabbed(os, "mass");
    writeTabbed(os, "Dmax");
    writeTabbed(os, "D10");
    writeTabbed(os, "D32");
    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cloudInfo::cloudInfo
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjectFile(obr, name),
    name_(name),
    obr_(obr),
    active_(true),
    log_(true),
    cloudNames_(),
    filePtrs_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cloudInfo::~cloudInfo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cloudInfo::read(const dictionary& dict)
{
    if (active_)
    {
        functionObjectFile::read(dict);

        log_ = dict.lookupOrDefault<Switch>("log", true);
        dict.lookup("clouds") >> cloudNames_;

        if (log_)
        {
            Info<< type() << " " << name_ << ": ";

            if (cloudNames_.size())
            {
                Info<< "applying to clouds:" << nl;
                forAll(cloudNames_, i)
                {
                    Info<< "    " << cloudNames_[i] << nl;
                }
                Info<< endl;
            }
            else
            {
                Info<< "no clouds to be processed" << nl << endl;
            }
        }

        if (writeToFile())
        {
            filePtrs_.setSize(cloudNames_.size());
            filePtrs_.clear();
            forAll(filePtrs_, fileI)
            {
                const word& cloudName = cloudNames_[fileI];
                filePtrs_.set(fileI, createFile(cloudName));
                writeFileHeader(filePtrs_[fileI]);
            }
        }
    }
}


void Foam::cloudInfo::execute()
{
    // Do nothing
}


void Foam::cloudInfo::end()
{
    // Do nothing
}


void Foam::cloudInfo::timeSet()
{
    // Do nothing
}


void Foam::cloudInfo::write()
{
    if (active_)
    {
        forAll(cloudNames_, cloudI)
        {
            const word& cloudName = cloudNames_[cloudI];

            const kinematicCloud& cloud =
                obr_.lookupObject<kinematicCloud>(cloudName);

            label nParcels = returnReduce(cloud.nParcels(), sumOp<label>());
            scalar massInSystem =
                returnReduce(cloud.massInSystem(), sumOp<scalar>());

            scalar Dmax = cloud.Dmax();
            scalar D10 = cloud.Dij(1, 0);
            scalar D32 = cloud.Dij(3, 2);

            if (Pstream::master())
            {
                writeTime(filePtrs_[cloudI]);
                filePtrs_[cloudI]
                    << token::TAB
                    << nParcels << token::TAB
                    << massInSystem << token::TAB
                    << Dmax << token::TAB
                    << D10 << token::TAB
                    << D32 << token::TAB
                    << endl;
            }

            if (log_) Info
                << type() << " " << name_ <<  " output:" << nl
                << "    number of parcels : " << nParcels << nl
                << "    mass in system    : " << massInSystem << nl
                << "    maximum diameter  : " << Dmax << nl
                << "    D10 diameter      : " << D10 << nl
                << "    D32 diameter      : " << D32 << nl
                << endl;
        }
    }
}


// ************************************************************************* //
