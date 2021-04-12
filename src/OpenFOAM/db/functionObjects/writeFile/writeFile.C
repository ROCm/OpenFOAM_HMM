/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2018 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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
#include "functionObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::label Foam::functionObjects::writeFile::addChars = 8;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::writeFile::initStream(Ostream& os) const
{
    os.setf(ios_base::scientific, ios_base::floatfield);
    os.precision(writePrecision_);
    os.width(charWidth());
}


Foam::fileName Foam::functionObjects::writeFile::baseFileDir() const
{
    // Put in undecomposed case
    // (Note: gives problems for distributed data running)

    fileName baseDir
    (
        fileObr_.time().globalPath()
      / functionObject::outputPrefix
    );

    // Append mesh name if not default region
    if (isA<polyMesh>(fileObr_))
    {
        const polyMesh& mesh = refCast<const polyMesh>(fileObr_);
        if (mesh.name() != polyMesh::defaultRegion)
        {
            baseDir /= mesh.name();
        }
    }
    baseDir.clean();  // Remove unneeded ".."

    return baseDir;
}


Foam::fileName Foam::functionObjects::writeFile::baseTimeDir() const
{
    return baseFileDir()/prefix_/fileObr_.time().timeName();
}


Foam::autoPtr<Foam::OFstream> Foam::functionObjects::writeFile::createFile
(
    const word& name,
    scalar timeValue
) const
{
    autoPtr<OFstream> osPtr;

    if (Pstream::master() && writeToFile_)
    {
        if (useUserTime_)
        {
            timeValue = fileObr_.time().timeToUserTime(timeValue);
        }

        const word timeName = Time::timeName(timeValue);

        fileName outputDir(baseFileDir()/prefix_/timeName);

        mkDir(outputDir);

        word fName(name);

        // Check if file already exists
        IFstream is(outputDir/(fName + ".dat"));
        if (is.good())
        {
            fName = fName + "_" + timeName;
        }

        osPtr.reset(new OFstream(outputDir/(fName + ".dat")));

        if (!osPtr->good())
        {
            FatalIOErrorInFunction(osPtr()) << "Cannot open file"
                << exit(FatalIOError);
        }

        initStream(osPtr());
    }

    return osPtr;
}


Foam::autoPtr<Foam::OFstream> Foam::functionObjects::writeFile::createFile
(
    const word& name
) const
{
    return createFile(name, startTime_);

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

Foam::functionObjects::writeFile::writeFile(const writeFile& wf)
:
    fileObr_(wf.fileObr_),
    prefix_(wf.prefix_),
    fileName_(wf.fileName_),
    filePtr_(nullptr),
    writePrecision_(wf.writePrecision_),
    writeToFile_(wf.writeToFile_),
    updateHeader_(wf.updateHeader_),
    writtenHeader_(wf.writtenHeader_),
    useUserTime_(wf.useUserTime_),
    startTime_(wf.startTime_)
{}


Foam::functionObjects::writeFile::writeFile
(
    const objectRegistry& obr,
    const fileName& prefix,
    const word& name,
    const bool writeToFile
)
:
    fileObr_(obr),
    prefix_(prefix),
    fileName_(name),
    filePtr_(nullptr),
    writePrecision_(IOstream::defaultPrecision()),
    writeToFile_(writeToFile),
    updateHeader_(true),
    writtenHeader_(false),
    useUserTime_(true),
    startTime_(obr.time().startTime().value())
{}


Foam::functionObjects::writeFile::writeFile
(
    const objectRegistry& obr,
    const fileName& prefix,
    const word& name,
    const dictionary& dict,
    const bool writeToFile
)
:
    writeFile(obr, prefix, name, writeToFile)
{
    read(dict);

    if (writeToFile_)
    {
        filePtr_ = createFile(fileName_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::writeFile::read(const dictionary& dict)
{
    writePrecision_ =
        dict.getCheckOrDefault
        (
            "writePrecision",
            IOstream::defaultPrecision(),
            labelMinMax::ge(0)
        );

    updateHeader_ = dict.getOrDefault("updateHeader", updateHeader_);

    // Only write on master
    writeToFile_ =
        Pstream::master() && dict.getOrDefault("writeToFile", writeToFile_);

    // Use user time, e.g. CA deg in preference to seconds
    useUserTime_ = dict.getOrDefault("useUserTime", true);

    return true;
}


Foam::OFstream& Foam::functionObjects::writeFile::file()
{
    if (!writeToFile_)
    {
        return Snull;
    }

    if (!filePtr_)
    {
        FatalErrorInFunction
            << "File pointer not allocated\n";
    }

    return *filePtr_;
}


bool Foam::functionObjects::writeFile::writeToFile() const
{
    return writeToFile_;
}


bool Foam::functionObjects::writeFile::canWriteHeader() const
{
    return writeToFile_ && (updateHeader_ || !writtenHeader_);
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


void Foam::functionObjects::writeFile::writeCurrentTime(Ostream& os) const
{
    const scalar timeValue =
    (
        useUserTime_
      ? fileObr_.time().timeOutputValue()
      : fileObr_.time().value()
    );

    os  << setw(charWidth()) << Time::timeName(timeValue);
}


void Foam::functionObjects::writeFile::writeBreak(Ostream& os) const
{
    writeHeader(os, "===");
}


// ************************************************************************* //
