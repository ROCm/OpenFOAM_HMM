/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2012-2018 OpenFOAM Foundation
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

    fileName baseDir =
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
    const scalar time0
) const
{
    autoPtr<OFstream> osPtr;

    if (Pstream::master() && writeToFile_)
    {
        const scalar time = useUserTime_ ?
            fileObr_.time().timeToUserTime(time0)
          : time0;

        const word timeName = Time::timeName(time);

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
    writtenHeader_(false),
    useUserTime_(true),
    startTime_(obr.time().startTime().value())
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
    writtenHeader_(false),
    useUserTime_(true),
    startTime_(obr.time().startTime().value())
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
        dict.lookupOrDefault("writePrecision", IOstream::defaultPrecision());

    // Only write on master process
    writeToFile_ = dict.lookupOrDefault("writeToFile", true);
    writeToFile_ = writeToFile_ && Pstream::master();

    // Use user time, e.g. CA deg in preference to seconds
    useUserTime_ = dict.lookupOrDefault("useUserTime", true);

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

    return *filePtr_;
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
    scalar timeNow = useUserTime_ ?
        fileObr_.time().timeOutputValue()
      : fileObr_.time().value();

    os  << setw(charWidth()) << Time::timeName(timeNow);
}


void Foam::functionObjects::writeFile::writeBreak(Ostream& os) const
{
    writeHeader(os, "===");
}


// ************************************************************************* //
