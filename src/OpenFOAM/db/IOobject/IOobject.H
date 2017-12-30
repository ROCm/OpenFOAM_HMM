/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

Class
    Foam::IOobject

Description
    IOobject defines the attributes of an object for which implicit
    objectRegistry management is supported, and provides the infrastructure
    for performing stream I/O.

    An IOobject is constructed with an object name, a class name, an instance
    path, a reference to a objectRegistry, and parameters determining its
    storage status.

    \par Read options

    Define what is done on object construction and explicit reads:
      - \par MUST_READ
        Object must be read from Istream on construction. \n
        Error if Istream does not exist or can't be read.
        Does not check timestamp or re-read.
      - \par MUST_READ_IF_MODIFIED
        Object must be read from Istream on construction. \n
        Error if Istream does not exist or can't be read. If object is
        registered its timestamp will be checked every timestep and possibly
        re-read.
      - \par READ_IF_PRESENT
        Read object from Istream if Istream exists, otherwise don't. \n
        Error only if Istream exists but can't be read.
        Does not check timestamp or re-read.
      - \par NO_READ
        Don't read

    \par Write options

    Define what is done on object destruction and explicit writes:
      - \par AUTO_WRITE
        Object is written automatically when requested to by the
        objectRegistry.
      - \par NO_WRITE
        No automatic write on destruction but can be written explicitly

SourceFiles
    IOobject.C
    IOobjectReadHeader.C
    IOobjectWriteHeader.C
    IOobjectPrint.C

\*---------------------------------------------------------------------------*/

#ifndef IOobject_H
#define IOobject_H

#include "fileName.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "InfoProxy.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Time;
class objectRegistry;

/*---------------------------------------------------------------------------*\
                           Class IOobject Declaration
\*---------------------------------------------------------------------------*/

class IOobject
{

public:

    // Public data types

        //- Enumeration defining the valid states of an IOobject
        enum objectState
        {
            GOOD,
            BAD
        };

        //- Enumeration defining the read options
        enum readOption
        {
            MUST_READ,
            MUST_READ_IF_MODIFIED,
            READ_IF_PRESENT,
            NO_READ
        };

        //- Enumeration defining the write options
        enum writeOption
        {
            AUTO_WRITE = 0,
            NO_WRITE = 1
        };

        //- Enumeration defining the file checking options
        enum fileCheckTypes
        {
            timeStamp,
            timeStampMaster,
            inotify,
            inotifyMaster
        };

        static const Enum<fileCheckTypes> fileCheckTypesNames;


private:

    // Private data

        //- Name
        word name_;

        //- Class name read from header
        word headerClassName_;

        //- Optional note
        string note_;

        //- Instance path component
        fileName instance_;

        //- Local path component
        fileName local_;

        //- objectRegistry reference
        const objectRegistry& db_;

        //- Read option
        readOption rOpt_;

        //- Write option
        writeOption wOpt_;

        //- Register object created from this IOobject with registry if true
        bool registerObject_;

        //- Is object same for all processors
        bool globalObject_;

        //- IOobject state
        objectState objState_;


protected:

    // Protected Member Functions

        //- Set the object state to bad
        void setBad(const string& s);


public:

    //- Runtime type information
    TypeName("IOobject");


    // Static data members

        //- Type of file modification checking
        static fileCheckTypes fileModificationChecking;


    // Static Member Functions

        //- Split path into instance, local, name components
        //
        // The splitting behaviour is as follows:
        // \verbatim
        //    input       | instance   | local | name
        //    ----------- | ---------- | ----- | ----
        //    a           |            |       | a
        //    a/b         | a          |       | b
        //    a/b/c/d     | a          | b/c   | d
        //    /a/b/c      | /a/b       |       | c
        //    ./a/b/c     | PWD/a/b    |       | c
        //    ../a/b/c    | PWD/../a/b |       | c
        //    a/b/        | ERROR      |       |
        // \endverbatim
        // where PWD is the Foam::cwd() current working directory
        static bool fileNameComponents
        (
            const fileName& path,
            fileName& instance,
            fileName& local,
            word& name
        );

        //- Create dot-delimited name.group
        template<class StringType>
        static inline word groupName(StringType name, const word& group);


    // Constructors

        //- Construct from name, instance, registry, io options
        IOobject
        (
            const word& name,
            const fileName& instance,
            const objectRegistry& registry,
            readOption r=NO_READ,
            writeOption w=NO_WRITE,
            bool registerObject=true
        );

        //- Construct from name, instance, local, registry, io options
        IOobject
        (
            const word& name,
            const fileName& instance,
            const fileName& local,
            const objectRegistry& registry,
            readOption r=NO_READ,
            writeOption w=NO_WRITE,
            bool registerObject=true,
            bool globalObject = false
        );

        //- Construct from path, registry, io options.
        //  Uses fileNameComponents() to split path into components.
        //  A path that starts with a '/' is regarded as a file system path.
        //  Paths starting with either './' or '../' are relative to
        //  current working directory (and replaced with absolute equivalents).
        //  All other paths are considered to be relative to the case.
        IOobject
        (
            const fileName& path,
            const objectRegistry& registry,
            readOption r=NO_READ,
            writeOption w=NO_WRITE,
            bool registerObject=true,
            bool globalObject = false
        );

        //- Construct as copy resetting registry
        IOobject
        (
            const IOobject& io,
            const objectRegistry& registry
        );

        //- Construct as copy resetting name
        IOobject
        (
            const IOobject& io,
            const word& name
        );

        //- Clone
        autoPtr<IOobject> clone() const
        {
            return autoPtr<IOobject>(new IOobject(*this));
        }

        //- Clone resetting registry
        autoPtr<IOobject> clone(const objectRegistry& registry) const
        {
            return autoPtr<IOobject>(new IOobject(*this, registry));
        }


    //- Destructor
    virtual ~IOobject();


    // Member Functions

      // General access

        //- Return the local objectRegistry
        const objectRegistry& db() const;

        //- Return time
        const Time& time() const;

        //- Return name
        inline const word& name() const;

        //- Return name of the class name read from header
        inline const word& headerClassName() const;

        //- Return non-constant access to the class name read from header
        inline word& headerClassName();

        //- Return the optional note
        inline const string& note() const;

        //- Return non-constant access to the optional note
        inline string& note();

        //- Rename
        virtual void rename(const word& newName)
        {
            name_ = newName;
        }

        //- Register object created from this IOobject with registry if true
        inline bool registerObject() const;

        //- Register object created from this IOobject with registry if true
        inline bool& registerObject();

        //- Is object same for all processors
        inline bool globalObject() const;

        //- Is object same for all processors
        inline bool& globalObject();


      // Read/write options

        //- The read option
        inline readOption readOpt() const;

        //- Non-constant access to the read option
        inline readOption& readOpt();

        //- The write option
        inline writeOption writeOpt() const;

        //- Non-constant access to the write option
        inline writeOption& writeOpt();


     // Path components

        //- Return group (extension part of name)
        word group() const;

        //- Return member (name without the extension)
        word member() const;

        const fileName& rootPath() const;

        const fileName& caseName() const;

        inline const fileName& instance() const;

        inline fileName& instance();

        inline const fileName& local() const;

        //- Return complete path
        fileName path() const;

        //- Return complete path with alternative instance and local
        fileName path
        (
            const word& instance,
            const fileName& local = fileName::null
        ) const;

        //- Return complete path + object name
        inline fileName objectPath() const;

        //- Helper for filePath that searches locally.
        //  When search is false, simply use the current instance,
        //  otherwise search previous instances.
        fileName localFilePath
        (
            const word& typeName,
            const bool search=true
        ) const;

        //- Helper for filePath that searches up if in parallel
        //  When search is false, simply use the current instance,
        //  otherwise search previous instances.
        fileName globalFilePath
        (
            const word& typeName,
            const bool search=true
        ) const;


      // Reading

        //- Read header
        bool readHeader(Istream& is);

        //- Read header (uses typeFilePath to find file) and check its info.
        //  Optionally checks headerClassName against the type-name.
        //  When search is false, simply use the current instance,
        //  otherwise search previous instances.
        template<class Type>
        bool typeHeaderOk
        (
            const bool checkType = true,
            const bool search = true,
            const bool verbose = true
        );

        //- Helper: warn that type does not support re-reading
        template<class Type>
        void warnNoRereading() const;


      // Writing

        //- Write the standard OpenFOAM file/dictionary banner
        //  Optionally without -*- C++ -*- editor hint (eg, for logs)
        static Ostream& writeBanner(Ostream& os, bool noHint=false);

        //- Write the standard file section divider
        static Ostream& writeDivider(Ostream& os);

        //- Write the standard end file divider
        static Ostream& writeEndDivider(Ostream& os);

        //- Write header
        bool writeHeader(Ostream& os) const;

        //- Write header. Allow override of type
        bool writeHeader(Ostream& os, const word& objectType) const;


      // Error Handling

        inline bool good() const;

        inline bool bad() const;


      // Info

        //- Return info proxy.
        //  Used to print token information to a stream
        InfoProxy<IOobject> info() const
        {
            return *this;
        }


    // Member operators

        void operator=(const IOobject& io);
};


template<>
Ostream& operator<<(Ostream& os, const InfoProxy<IOobject>& ip);


//- Template function for obtaining global vs. local status
template<class T>
inline bool typeGlobal()
{
    return false;
}


//- Template function for obtaining local or global filePath
template<class T>
inline fileName typeFilePath(const IOobject& io, const bool search = true)
{
    return
    (
        typeGlobal<T>()
      ? io.globalFilePath(T::typeName, search)
      : io.localFilePath(T::typeName, search)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IOobjectI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "IOobjectTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
