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

Class
    Foam::IOobject

Description
    Defines the attributes of an object for which implicit
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

    When serializing, the IOobject characteristics are typically written
    as a \c FoamFile header, which is a sub-dictionary with the following
    type of content:

    \table
        Property | Description                            | Type | Reqd | Deflt
        version  | The base format version                | float | no  | 2.0
        format   | The stream format (ascii/binary)       | word  | yes |
        arch     | The architecture string                | string | no |
        note     | Descriptive note about the object      | string | no |
        location | The relative location of the object    | string | no |
        class    | The type of the object                 | word | yes |
        object   | The name of the object                 | word | yes |
    \endtable

Note
    Specifying registered does not result in the IOobject itself being
    registered. It is only serves as guidance for a regIOobject using it.

See also
    Foam::objectRegistry
    Foam::regIOobject

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
#include "IOstreamOption.H"
#include "InfoProxy.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Time;
class dictionary;
class objectRegistry;

/*---------------------------------------------------------------------------*\
                           Class IOobject Declaration
\*---------------------------------------------------------------------------*/

class IOobject
{
public:

    // Public Data Types

        //- Enumeration defining the valid states of an IOobject
        enum objectState : char
        {
            GOOD,
            BAD
        };

        //- Enumeration defining the read options
        enum readOption : char
        {
            MUST_READ,
            MUST_READ_IF_MODIFIED,
            READ_IF_PRESENT,
            NO_READ
        };

        //- Enumeration defining the write options
        enum writeOption : char
        {
            AUTO_WRITE = 0,
            NO_WRITE = 1
        };

        //- Enumeration defining the file checking options
        enum fileCheckTypes : char
        {
            timeStamp,
            timeStampMaster,
            inotify,
            inotifyMaster
        };

        //- Names for the fileCheckTypes
        static const Enum<fileCheckTypes> fileCheckTypesNames;


private:

    // Static Data Members

        //- Use an output file banner, enabled by default
        static bool bannerEnabled_;


    // Private Data

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

        //- Read option
        readOption rOpt_;

        //- Write option
        writeOption wOpt_;

        //- Should object created with this IOobject be registered?
        bool registerObject_;

        //- Is object same for all processors?
        bool globalObject_;

        //- IOobject state
        objectState objState_;

        //- The sizeof (label) in bytes, possibly read from the header
        unsigned char sizeofLabel_;

        //- The sizeof (scalar) in bytes, possibly read from the header
        unsigned char sizeofScalar_;

        //- Reference to the objectRegistry
        const objectRegistry& db_;


protected:

    // Protected Member Functions

        //- Helper: write content for FoamFile IOobject header
        //- with optional meta information.
        static void writeHeaderContent
        (
            Ostream& os,
            const IOobject& io,
            const word& objectType,
            const dictionary* metaDataDict = nullptr
        );

        //- Helper: write dictionary content for FoamFile header
        //- with optional meta information.
        static void writeHeaderContent
        (
            dictionary& dict,
            const IOobject& io,
            const word& objectType,
            IOstreamOption streamOpt,
            const dictionary* metaDataDict = nullptr
        );

        //- Set the object state to bad
        void setBad(const string& s);


public:

    //- Declare type-name, virtual type (with debug switch)
    TypeName("IOobject");


    // Static Data

        //- Character for scoping object names (':' or '_')
        //  Change with caution.
        static char scopeSeparator;

        //- Type of file modification checking
        static fileCheckTypes fileModificationChecking;

        //- Time skew (seconds) for file modification checks
        static float fileModificationSkew;

        //- Max number of times to poll for file modification changes
        static int maxFileModificationPolls;


    // Static Functions

        //- Status of output file banner
        static bool bannerEnabled() noexcept
        {
            return bannerEnabled_;
        }

        //- Enable/disable an output file banner
        //  \return the previous value
        static bool bannerEnabled(bool on) noexcept
        {
            bool old(bannerEnabled_);
            bannerEnabled_ = on;
            return old;
        }

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

        //- Create dot-delimited name.group string
        //  An empty group is ignored.
        template<class StringType>
        static inline word groupName(StringType base, const word& group);

        //- Return group (extension part of name)
        static word group(const word& name);

        //- Return member (name without the extension)
        static word member(const word& name);

        //- Create scope:name or scope_name string
        //  An empty scope is ignored.
        static inline word scopedName
        (
            const std::string& scope,
            const word& name
        );

        //- Return the IOobject, but also consider an alternative file name.
        //
        //  \param io The expected IOobject to use
        //  \param altFile Alternative fileName (ignored if empty).
        //  \param ioName  The alternative name for the IOobject when
        //      the altFile resolves to a directory.
        //
        //  \note If the alternative fileName is a non-empty string,
        //      it defines the location but uses all other properties of the
        //      expected IOobject.
        //      The location may be an absolute or a relative path.
        //      If it corresponds to a directory, the name of the
        //      expected IOobject will be used in its resolution.
        //      This expected name can provided via the ioName parameter.
        static IOobject selectIO
        (
            const IOobject& io,
            const fileName& altFile,
            const word& ioName = ""
        );


    // Generated Methods

        //- Copy construct
        IOobject(const IOobject&) = default;

        //- Destructor
        virtual ~IOobject() = default;


    // Constructors

        //- Construct from name, instance, registry, io options
        IOobject
        (
            const word& name,
            const fileName& instance,
            const objectRegistry& registry,
            readOption r=NO_READ,
            writeOption w=NO_WRITE,
            bool registerObject = true,
            bool globalObject = false
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
            bool registerObject = true,
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
            bool registerObject = true,
            bool globalObject = false
        );

        //- Copy construct, resetting registry
        IOobject
        (
            const IOobject& io,
            const objectRegistry& registry
        );

        //- Copy construct, resetting name
        IOobject
        (
            const IOobject& io,
            const word& name
        );

        //- Copy construct, resetting io options
        IOobject
        (
            const IOobject& io,
            readOption,
            writeOption
        );

        //- Clone
        autoPtr<IOobject> clone() const
        {
            return autoPtr<IOobject>::New(*this);
        }

        //- Clone resetting registry
        autoPtr<IOobject> clone(const objectRegistry& registry) const
        {
            return autoPtr<IOobject>::New(*this, registry);
        }


    // Member Functions

    // General Access

        //- Return the local objectRegistry
        const objectRegistry& db() const noexcept;

        //- Return Time associated with the objectRegistry
        const Time& time() const;

        //- Return name
        inline const word& name() const noexcept;

        //- Return name of the class name read from header
        inline const word& headerClassName() const noexcept;

        //- Return non-constant access to the class name read from header
        inline word& headerClassName() noexcept;

        //- Return the optional note
        inline const string& note() const noexcept;

        //- Return non-constant access to the optional note
        inline string& note() noexcept;

        //- Rename
        virtual void rename(const word& newName)
        {
            name_ = newName;
        }

        //- Should object created with this IOobject be registered?
        inline bool registerObject() const noexcept;

        //- Change registration preference, return previous value
        inline bool registerObject(bool on) noexcept;

        //- Is object same for all processors?
        inline bool globalObject() const noexcept;

        //- Change global-object status, return previous value
        inline bool globalObject(bool on) noexcept;

        //- The sizeof (label) in bytes, possibly read from the header
        inline unsigned labelByteSize() const noexcept;

        //- The sizeof (scalar) in bytes, possibly read from the header
        inline unsigned scalarByteSize() const noexcept;


    // Checks

        //- Test if headerClassName() equals the given class name
        inline bool isHeaderClassName(const word& clsName) const;

        //- Test if headerClassName() equals Type::typeName
        template<class Type>
        inline bool isHeaderClassName() const;


    // Meta-data

        //- Return pointer to meta-data (if any) or nullptr
        virtual const dictionary* findMetaData() const noexcept;


    // Read/write options

        //- The read option
        inline readOption readOpt() const noexcept;

        //- Change the read option, return previous value
        inline readOption readOpt(readOption opt) noexcept;

        //- The write option
        inline writeOption writeOpt() const noexcept;

        //- Change the write option, return previous value
        inline writeOption writeOpt(writeOption opt) noexcept;


    // Path components

        //- Return group (extension part of name)
        inline word group() const;

        //- Return member (name without the extension)
        inline word member() const;

        const fileName& rootPath() const;

        const fileName& caseName() const;

        inline const fileName& instance() const noexcept;

        inline fileName& instance() noexcept;

        inline const fileName& local() const noexcept;

        //- The complete path
        fileName path() const;

        //- The complete path with alternative instance and local
        fileName path
        (
            const word& instance,
            const fileName& local = fileName::null
        ) const;

        //- The complete path + object name
        inline fileName objectPath() const;

        //- The object path relative to the root
        fileName objectRelPath() const;

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

        //- Parse 'FoamFile' header contents and set the IOobject
        //- characteristics and return the stream characteristics.
        IOstreamOption parseHeader(const dictionary& headerDict);

        //- Read header ('FoamFile' dictionary) and set the
        //- IOobject and stream characteristics.
        bool readHeader(Istream& is);

        //- Read header (the 'FoamFile' dictionary) and set the
        //- IOobject and stream characteristics.
        //  Saves the header content in the given dictionary.
        bool readHeader(dictionary& headerDict, Istream& is);

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
        //  Optionally without editor syntax hint (eg, for logs)
        static Ostream& writeBanner(Ostream& os, const bool noSyntaxHint=false);

        //- Write the standard file section divider
        static Ostream& writeDivider(Ostream& os);

        //- Write the standard end file divider
        static Ostream& writeEndDivider(Ostream& os);

        //- Write header with current type()
        bool writeHeader(Ostream& os) const;

        //- Write header with override of type
        bool writeHeader(Ostream& os, const word& objectType) const;

        //- Write header into a dictionary with current type()
        //- and given output format
        void writeHeader(dictionary& dict, IOstreamOption streamOpt) const;

        //- Write header into a dictionary with override of type
        //- and given output format
        void writeHeader
        (
            dictionary& dict,
            const word& objectType,
            IOstreamOption streamOpt
        ) const;


    // Error Handling

        inline bool good() const noexcept;

        inline bool bad() const noexcept;


    // Info

        //- Return info proxy.
        //  Used to print token information to a stream
        InfoProxy<IOobject> info() const
        {
            return *this;
        }


    // Member Operators

        void operator=(const IOobject& io);


    // Housekeeping

        //- Access to the read option
        //  \deprecated(2021-03) - use readOpt(readOption)
        readOption& readOpt() noexcept { return rOpt_; }

        //- Access to the write option
        //  \deprecated(2021-03) - use writeOpt(writeOption)
        writeOption& writeOpt() noexcept { return wOpt_; }

        //- Access to the register object option
        //  \deprecated(2021-03) - use registerObject(bool)
        bool& registerObject() noexcept { return registerObject_; }

        //- Access to the global object option
        //  \deprecated(2021-03) - use globalObject(bool)
        bool& globalObject() noexcept { return globalObject_; }
};


//- Specialization for \c void always returns true (no headerClassName check).
template<>
inline bool IOobject::isHeaderClassName<void>() const
{
    return true;
}


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
