/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "IOobject.H"
#include "Time.H"
#include "Istream.H"
#include "registerSwitch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IOobject, 0);
}

bool Foam::IOobject::bannerEnabled_(true);

char Foam::IOobject::scopeSeparator
(
    #ifdef _WIN32
    // Windows: using ':' causes scoping conflicts with d:/path etc
    Foam::debug::infoSwitch("scopeSeparator", '_')
    #else
    Foam::debug::infoSwitch("scopeSeparator", ':')
    #endif
);

const Foam::Enum
<
    Foam::IOobject::fileCheckTypes
>
Foam::IOobject::fileCheckTypesNames
({
    { fileCheckTypes::timeStamp, "timeStamp" },
    { fileCheckTypes::timeStampMaster, "timeStampMaster" },
    { fileCheckTypes::inotify, "inotify" },
    { fileCheckTypes::inotifyMaster, "inotifyMaster" },
});

// Default fileCheck type
Foam::IOobject::fileCheckTypes Foam::IOobject::fileModificationChecking
(
    fileCheckTypesNames.get
    (
        "fileModificationChecking",
        debug::optimisationSwitches()
    )
);


float Foam::IOobject::fileModificationSkew
(
    Foam::debug::floatOptimisationSwitch("fileModificationSkew", 30)
);
registerOptSwitch
(
    "fileModificationSkew",
    float,
    Foam::IOobject::fileModificationSkew
);

int Foam::IOobject::maxFileModificationPolls
(
    Foam::debug::optimisationSwitch("maxFileModificationPolls", 1)
);
registerOptSwitch
(
    "maxFileModificationPolls",
    int,
    Foam::IOobject::maxFileModificationPolls
);


//! \cond file-scope
namespace Foam
{
    // Register re-reader
    class addfileModificationCheckingToOpt
    :
        public ::Foam::simpleRegIOobject
    {
    public:

        addfileModificationCheckingToOpt
            (const addfileModificationCheckingToOpt&) = delete;

        void operator=
            (const addfileModificationCheckingToOpt&) = delete;

        explicit addfileModificationCheckingToOpt(const char* name)
        :
            ::Foam::simpleRegIOobject(Foam::debug::addOptimisationObject, name)
        {}

        virtual ~addfileModificationCheckingToOpt() = default;

        virtual void readData(Foam::Istream& is)
        {
            IOobject::fileModificationChecking =
                IOobject::fileCheckTypesNames.read(is);
        }

        virtual void writeData(Foam::Ostream& os) const
        {
            os <<  IOobject::fileCheckTypesNames
                [IOobject::fileModificationChecking];
        }
    };

    addfileModificationCheckingToOpt addfileModificationCheckingToOpt_
    (
        "fileModificationChecking"
    );

} // End namespace Foam
//! \endcond


// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

bool Foam::IOobject::fileNameComponents
(
    const fileName& path,
    fileName& instance,
    fileName& local,
    word& name
)
{
    // Convert explicit relative file-system path to absolute file-system path.
    if (path.starts_with("./") || path.starts_with("../"))
    {
        fileName absPath(cwd()/path);
        absPath.clean();  // Remove unneeded ".."

        return fileNameComponents(absPath, instance, local, name);
    }

    instance.clear();
    local.clear();
    name.clear();

    // Called with directory
    if (isDir(path))
    {
        WarningInFunction
            << " called with directory: " << path << endl;

        return false;
    }

    const auto first = path.find('/');
    const auto last  = path.rfind('/');

    // The raw length of name (without validating for word chars)
    auto nameLen = path.size();

    if (first == std::string::npos)
    {
        // No '/' found (or empty entirely)
        // => no instance or local

        name = word::validate(path);
    }
    else if
    (
        first == 0
        #ifdef _WIN32
     || (first == 2 && path[1] == ':')  // Eg, d:/path
        #endif
    )
    {
        // Absolute path (starts with '/' or 'd:/')
        // => no local

        instance = path.substr(0, last);

        const std::string ending = path.substr(last+1);
        nameLen = ending.size();  // The raw length of name
        name = word::validate(ending);
    }
    else
    {
        // Normal case.
        // First part is instance, remainder is local
        instance = path.substr(0, first);

        if (last > first)
        {
            // With local
            local = path.substr(first+1, last-first-1);
        }

        const std::string ending = path.substr(last+1);
        nameLen = ending.size();  // The raw length of name
        name = word::validate(ending);
    }

    // Check for valid (and stripped) name, regardless of the debug level
    if (!nameLen || nameLen != name.size())
    {
        WarningInFunction
            << "has invalid word for name: \"" << name
            << "\"\nwhile processing path: " << path << endl;

        return false;
    }

    return true;
}


Foam::IOobject Foam::IOobject::selectIO
(
    const IOobject& io,
    const fileName& altFile,
    const word& ioName
)
{
    if (altFile.empty())
    {
        return io;
    }

    // Construct from file path instead

    fileName altPath = altFile;

    if (isDir(altPath))
    {
        // Resolve directories as well

        if (ioName.empty())
        {
            altPath /= io.name();
        }
        else
        {
            altPath /= ioName;
        }
    }
    altPath.expand();


    return
        IOobject
        (
            altPath,
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            io.registerObject(),
            io.globalObject()
        );
}


Foam::word Foam::IOobject::group(const word& name)
{
    const auto i = name.rfind('.');

    if (i == std::string::npos || i == 0)
    {
        return word::null;
    }

    return name.substr(i+1);
}


Foam::word Foam::IOobject::member(const word& name)
{
    const auto i = name.rfind('.');

    if (i == std::string::npos || i == 0)
    {
        return name;
    }

    return name.substr(0, i);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IOobject::IOobject
(
    const word& name,
    const fileName& instance,
    const objectRegistry& registry,
    readOption ro,
    writeOption wo,
    bool registerObject,
    bool globalObject
)
:
    name_(name),
    headerClassName_(typeName),
    note_(),
    instance_(instance),
    local_(),
    rOpt_(ro),
    wOpt_(wo),
    registerObject_(registerObject),
    globalObject_(globalObject),
    objState_(GOOD),
    sizeofLabel_(static_cast<unsigned char>(sizeof(label))),
    sizeofScalar_(static_cast<unsigned char>(sizeof(scalar))),

    db_(registry)
{
    if (objectRegistry::debug)
    {
        InfoInFunction
            << "Constructing IOobject called " << name_
            << " of type " << headerClassName_
            << endl;
    }
}


Foam::IOobject::IOobject
(
    const word& name,
    const fileName& instance,
    const fileName& local,
    const objectRegistry& registry,
    readOption ro,
    writeOption wo,
    bool registerObject,
    bool globalObject
)
:
    name_(name),
    headerClassName_(typeName),
    note_(),
    instance_(instance),
    local_(local),
    rOpt_(ro),
    wOpt_(wo),
    registerObject_(registerObject),
    globalObject_(globalObject),
    objState_(GOOD),
    sizeofLabel_(static_cast<unsigned char>(sizeof(label))),
    sizeofScalar_(static_cast<unsigned char>(sizeof(scalar))),

    db_(registry)
{
    if (objectRegistry::debug)
    {
        InfoInFunction
            << "Constructing IOobject called " << name_
            << " of type " << headerClassName_
            << endl;
    }
}


Foam::IOobject::IOobject
(
    const fileName& path,
    const objectRegistry& registry,
    readOption ro,
    writeOption wo,
    bool registerObject,
    bool globalObject
)
:
    name_(),
    headerClassName_(typeName),
    note_(),
    instance_(),
    local_(),
    rOpt_(ro),
    wOpt_(wo),
    registerObject_(registerObject),
    globalObject_(globalObject),
    objState_(GOOD),
    sizeofLabel_(static_cast<unsigned char>(sizeof(label))),
    sizeofScalar_(static_cast<unsigned char>(sizeof(scalar))),

    db_(registry)
{
    if (!fileNameComponents(path, instance_, local_, name_))
    {
        FatalErrorInFunction
            << " invalid path specification"
            << exit(FatalError);
    }

    if (objectRegistry::debug)
    {
        InfoInFunction
            << "Constructing IOobject called " << name_
            << " of type " << headerClassName_
            << endl;
    }
}


Foam::IOobject::IOobject
(
    const IOobject& io,
    const objectRegistry& registry
)
:
    name_(io.name_),
    headerClassName_(io.headerClassName_),
    note_(io.note_),
    instance_(io.instance_),
    local_(io.local_),
    rOpt_(io.rOpt_),
    wOpt_(io.wOpt_),
    registerObject_(io.registerObject_),
    globalObject_(io.globalObject_),
    objState_(io.objState_),
    sizeofLabel_(io.sizeofLabel_),
    sizeofScalar_(io.sizeofScalar_),

    db_(registry)
{}


Foam::IOobject::IOobject
(
    const IOobject& io,
    const word& name
)
:
    name_(name),
    headerClassName_(io.headerClassName_),
    note_(io.note_),
    instance_(io.instance_),
    local_(io.local_),
    rOpt_(io.rOpt_),
    wOpt_(io.wOpt_),
    registerObject_(io.registerObject_),
    globalObject_(io.globalObject_),
    objState_(io.objState_),
    sizeofLabel_(io.sizeofLabel_),
    sizeofScalar_(io.sizeofScalar_),

    db_(io.db_)
{}


Foam::IOobject::IOobject
(
    const IOobject& io,
    readOption ro,
    writeOption wo
)
:
    IOobject(io)
{
    rOpt_ = ro;
    wOpt_ = wo;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::objectRegistry& Foam::IOobject::db() const noexcept
{
    return db_;
}


const Foam::Time& Foam::IOobject::time() const
{
    return db_.time();
}


const Foam::fileName& Foam::IOobject::rootPath() const
{
    return time().rootPath();
}


const Foam::fileName& Foam::IOobject::caseName() const
{
    return time().caseName();
}


Foam::fileName Foam::IOobject::path() const
{
    // A file is 'outside' of the case if it has been specified using an
    // absolute path

    const auto first = instance().find('/');

    if
    (
        first == 0
        #ifdef _WIN32
     || (first == 2 && instance()[1] == ':')  // Eg, d:/path
        #endif
    )
    {
        // Absolute path (starts with '/' or 'd:/')
        return instance();
    }

    return rootPath()/caseName()/instance()/db_.dbDir()/local();
}


Foam::fileName Foam::IOobject::path
(
    const word& instance,
    const fileName& local
) const
{
    // Note: can only be called with relative instance since is word type
    return rootPath()/caseName()/instance/db_.dbDir()/local;
}


Foam::fileName Foam::IOobject::objectRelPath() const
{
    // A file is 'outside' of the case if it has been specified using an
    // absolute path

    const auto first = instance().find('/');

    if
    (
        first == 0
        #ifdef _WIN32
     || (first == 2 && instance()[1] == ':')  // Eg, d:/path
        #endif
    )
    {
        // Absolute path (starts with '/' or 'd:/')
        return instance()/name();
    }

    return instance()/db_.dbDir()/local()/name();
}


Foam::fileName Foam::IOobject::localFilePath
(
    const word& typeName,
    const bool search
) const
{
    // Do not check for undecomposed files
    return fileHandler().filePath(false, *this, typeName, search);
}


Foam::fileName Foam::IOobject::globalFilePath
(
    const word& typeName,
    const bool search
) const
{
    // Check for undecomposed files
    return fileHandler().filePath(true, *this, typeName, search);
}


void Foam::IOobject::setBad(const string& s)
{
    if (objState_ != GOOD)
    {
        FatalErrorInFunction
            << "Recurrent failure for object " << s
            << exit(FatalError);
    }

    if (error::level)
    {
        InfoInFunction
            << "Broken object " << s << info() << endl;
    }

    objState_ = BAD;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::IOobject::operator=(const IOobject& io)
{
    name_ = io.name_;
    headerClassName_ = io.headerClassName_;
    note_ = io.note_;
    instance_ = io.instance_;
    local_ = io.local_;
    rOpt_ = io.rOpt_;
    wOpt_ = io.wOpt_;
    globalObject_ = io.globalObject_;
    objState_ = io.objState_;
    sizeofLabel_ = io.sizeofLabel_;
    sizeofScalar_ = io.sizeofScalar_;
}


// ************************************************************************* //
