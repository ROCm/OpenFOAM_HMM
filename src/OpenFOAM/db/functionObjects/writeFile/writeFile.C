/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2017 OpenCFD Ltd.
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

#include "writeFile.H"
#include "Time.H"
#include "polyMesh.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::functionObjects::writeFile::outputPrefix
(
    "postProcessing"
);

Foam::label Foam::functionObjects::writeFile::addChars = 7;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::writeFile::initStream(Ostream& os) const
{
    os.setf(ios_base::scientific, ios_base::floatfield);
    os.precision(writePrecision_);
    os.width(charWidth());
}


Foam::fileName Foam::functionObjects::writeFile::baseFileDir() const
{
    fileName baseDir = fileObr_.time().path();

    if (Pstream::parRun())
    {
        // Put in undecomposed case (Note: gives problems for
        // distributed data running)
        baseDir = baseDir/".."/outputPrefix;
    }
    else
    {
        baseDir = baseDir/outputPrefix;
    }

    // Append mesh name if not default region
    if (isA<polyMesh>(fileObr_))
    {
        const polyMesh& mesh = refCast<const polyMesh>(fileObr_);
        if (mesh.name() != polyMesh::defaultRegion)
        {
            baseDir = baseDir/mesh.name();
        }
    }

    // Remove any ".."
    baseDir.clean();

    return baseDir;
}


Foam::fileName Foam::functionObjects::writeFile::baseTimeDir() const
{
    return baseFileDir()/prefix_/fileObr_.time().timeName();
}


Foam::autoPtr<Foam::OFstream> Foam::functionObjects::writeFile::createFile
(
    const word& name
) const
{
    autoPtr<OFstream> osPtr;

    if (Pstream::master() && writeToFile_)
    {
        const scalar startTime = fileObr_.time().startTime().value();
        const scalar userStartTime = fileObr_.time().timeToUserTime(startTime);
        const word startTimeName = Time::timeName(userStartTime);

        fileName outputDir(baseFileDir()/prefix_/startTimeName);

        mkDir(outputDir);

        word fName(name);

        // Check if file already exists
        IFstream is(outputDir/(fName + ".dat"));
        if (is.good())
        {
            fName = fName + "_" + startTimeName;
        }

        osPtr.reset(new OFstream(outputDir/(fName + ".dat")));

        initStream(osPtr());
    }

    return osPtr;
}


void Foam::functionObjects::writeFile::resetFile(const word& fileName)
{
    fileName_ = fileName;
    filePtr_ = createFile(fileName_);
}


Foam::Omanip<int> Foam::functionObjects::writeFile::valueWidth
(
    const label offset
) const
{
    return setw(writePrecision_ + addChars + offset);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::writeFile::writeFile
(
    const objectRegistry& obr,
    const word& prefix
)
:
    fileObr_(obr),
    prefix_(prefix),
    fileName_("undefined"),
    filePtr_(),
    writePrecision_(IOstream::defaultPrecision()),
    writeToFile_(true),
    writtenHeader_(false)
{}


Foam::functionObjects::writeFile::writeFile
(
    const objectRegistry& obr,
    const word& prefix,
    const word& fileName,
    const dictionary& dict
)
:
    fileObr_(obr),
    prefix_(prefix),
    fileName_(fileName),
    filePtr_(),
    writePrecision_(IOstream::defaultPrecision()),
    writeToFile_(true),
    writtenHeader_(false)
{
    read(dict);

    if (writeToFile_)
    {
        filePtr_ = createFile(fileName_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::writeFile::~writeFile()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::writeFile::read(const dictionary& dict)
{
    writePrecision_ =
        dict.lookupOrDefault("writePrecision", IOstream::defaultPrecision());

    // Only write on master process
    writeToFile_ = dict.lookupOrDefault("writeToFile", true);
    writeToFile_ = writeToFile_ && Pstream::master();

    return true;
}


Foam::OFstream& Foam::functionObjects::writeFile::file()
{
    if (!writeToFile_)
    {
        return Snull;
    }

    if (!filePtr_.valid())
    {
        FatalErrorInFunction
            << "File pointer not allocated";
    }

    return filePtr_();
}


bool Foam::functionObjects::writeFile::writeToFile() const
{
    return writeToFile_;
}


Foam::label Foam::functionObjects::writeFile::charWidth() const
{
    return writePrecision_ + addChars;
}


void Foam::functionObjects::writeFile::writeCommented
(
    Ostream& os,
    const string& str
) const
{
    os  << setw(1) << "#";

    if (str.size())
    {
        os  << setw(1) << ' '
            << setf(ios_base::left) << setw(charWidth() - 2) << str.c_str();
    }
}


void Foam::functionObjects::writeFile::writeTabbed
(
    Ostream& os,
    const string& str
) const
{
    os  << tab << setw(charWidth()) << str.c_str();
}


void Foam::functionObjects::writeFile::writeHeader
(
    Ostream& os,
    const string& str
) const
{
    writeCommented(os, str);
    os  << nl;
}


void Foam::functionObjects::writeFile::writeTime(Ostream& os) const
{
    const scalar timeNow = fileObr_.time().timeOutputValue();
    os  << setw(charWidth()) << Time::timeName(timeNow);
}


void Foam::functionObjects::writeFile::writeBreak(Ostream& os) const
{
    writeHeader(os, "===");
}


// ************************************************************************* //
