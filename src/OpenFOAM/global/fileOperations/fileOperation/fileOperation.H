/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::fileOperation

Description
    An encapsulation of filesystem-related operations.

Namespace
    Foam::fileOperations

Description
    Namespace for implementations of a fileOperation

\*---------------------------------------------------------------------------*/

#ifndef fileOperation_H
#define fileOperation_H

#include "ISstream.H"
#include "Ostream.H"
#include "autoPtr.H"
#include "fileNameList.H"
#include "instantList.H"
#include "fileMonitor.H"
#include "refPtr.H"
#include "Enum.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class IOobject;
class regIOobject;
class objectRegistry;
class Time;

/*---------------------------------------------------------------------------*\
                         Class fileOperation Declaration
\*---------------------------------------------------------------------------*/

class fileOperation
{
public:

        //- Enumeration for the location of an IOobject
        enum pathType : int
        {
            NOTFOUND = 0,           //!< Not found
            ABSOLUTE,               //!< instance is absolute directory
            OBJECT,                 //!< io.objectPath() exists
            WRITEOBJECT,            //!< write path exists

            // NOTE: increasing precedence (uncollated, collated, rank-collated)

            PROCUNCOLLATED,
                        //!< objectPath exists in 'processorN'
            PROCBASEOBJECT = PROCUNCOLLATED + 1,
                        //!< objectPath exists in 'processorsNN'
            PROCOBJECT = PROCBASEOBJECT + 1,
                        //!< objectPath exists in 'processorsNN_first-last'

            PARENTOBJECT,           //!< parent of object path
            FINDINSTANCE,           //!< file found in time directory
            PROCUNCOLLATEDINSTANCE, //!< as PROCUNCOLLATED but with instance
            PROCBASEINSTANCE,       //!< as PROCBASEOBJECT but with instance
            PROCINSTANCE            //!< as PROCOBJECT but with instance
        };
        static const Enum<pathType> pathTypeNames_;

        //- Augment fileName with pathType and local offset
        typedef Tuple2<fileName, Tuple2<pathType, int>> dirIndex;
        typedef List<dirIndex> dirIndexList;

        //- For addressing a range of processors,
        //- identical to UPstream::rangeType
        typedef IntRange<int> procRangeType;


protected:

    // Protected Data

        //- Communicator to use
        const label comm_;

        //- Distributed roots (parallel run)
        mutable bool distributed_;

        //- Detected processors directories
        mutable HashTable<dirIndexList> procsDirs_;

        //- File-change monitor for all registered files
        mutable autoPtr<fileMonitor> monitorPtr_;


   // Protected Member Functions

        //- Get or create fileMonitor singleton
        fileMonitor& monitor() const;

        //- Retrieve list of IO ranks from FOAM_IORANKS env variable
        static labelList ioRanks();

        //- Merge two times
        static void mergeTimes
        (
            const instantList& extraTimes,
            const word& constantName,
            instantList& times
        );

        //- Helper: check for file (isFile) or directory (!isFile)
        static bool isFileOrDir(const bool isFile, const fileName&);

        //- Lookup name of processorsDDD using cache.
        //  \return empty fileName if not found.
        refPtr<dirIndexList> lookupAndCacheProcessorsPath
        (
            const fileName& objectPath,
            const bool syncPar
        ) const;

        //- Lookup name of processorsDDD using cache.
        //  \note To be called on all processors
        //  \return empty fileName if not found.
        virtual refPtr<dirIndexList> lookupProcessorsPath
        (
            const fileName& objectPath
        ) const;

        //- Does IOObject exist.
        //- Is either a directory (empty name()) or a file
        bool exists(IOobject& io) const;


public:

    // Static Data

        //- Return the processors directory name (usually "processors")
        static word processorsBaseDir;

        //- Name of the default fileHandler
        static word defaultFileHandler;


    // Public Data Types

        //- Runtime type information
        TypeName("fileOperation");


    //- Static fileOperation
    static autoPtr<fileOperation> fileHandlerPtr_;

    //- Static construct the commonly used uncollatedFileOperation
    static autoPtr<fileOperation> NewUncollated();


    // Constructors

        //- Construct from communicator, optionally with distributed roots
        explicit fileOperation
        (
            const label comm,
            const bool distributedRoots = false
        );


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fileOperation,
            word,
            (
                bool verbose
            ),
            (verbose)
        );


    // Selectors

        //- Select fileHandler-type
        static autoPtr<fileOperation> New
        (
            const word& handlerType,
            bool verbose = false
        );


    //- Destructor
    virtual ~fileOperation() = default;


   // Static Functions

        //- Sort directory entries according to time value,
        //  with "constant" appearing first (if it exists)
        static instantList sortTimes
        (
            const fileNameList& dirEntries,
            const word& constantName = "constant"
        );


    // Member Functions

        //- Distributed roots (parallel run)
        bool distributed() const noexcept
        {
            return distributed_;
        }

        //- Set distributed roots on/off (mutable)
        //  \return old value
        bool distributed(bool on) const noexcept;


        // OSSpecific equivalents

            //- Make directory
            virtual bool mkDir(const fileName&, mode_t=0777) const = 0;

            //- Set the file mode
            virtual bool chMod(const fileName&, const mode_t) const = 0;

            //- Return the file mode
            virtual mode_t mode
            (
                const fileName&,
                const bool followLink = true
            ) const = 0;

            //- Return the file type: DIRECTORY, FILE or LINK
            virtual fileName::Type type
            (
                const fileName&,
                const bool followLink = true
            ) const = 0;

            //- Does the name exist (as DIRECTORY or FILE) in the file system?
            //  Optionally enable/disable check for gzip file.
            virtual bool exists
            (
                const fileName&,
                const bool checkGzip=true,
                const bool followLink = true
            ) const = 0;

            //- Does the name exist as a DIRECTORY in the file system?
            virtual bool isDir
            (
                const fileName&,
                const bool followLink = true
            ) const = 0;

            //- Does the name exist as a FILE in the file system?
            //  Optionally enable/disable check for gzip file.
            virtual bool isFile
            (
                const fileName&,
                const bool checkGzip=true,
                const bool followLink = true
            ) const = 0;

            //- Return size of file
            virtual off_t fileSize
            (
                const fileName&,
                const bool followLink = true
            ) const = 0;

            //- Return time of last file modification
            virtual time_t lastModified
            (
                const fileName&,
                const bool followLink = true
            ) const = 0;

            //- Return time of last file modification
            virtual double highResLastModified
            (
                const fileName&,
                const bool followLink = true
            ) const = 0;

            //- Read a directory and return the entries as a string list
            virtual fileNameList readDir
            (
                const fileName&,
                const fileName::Type=fileName::FILE,
                const bool filtergz=true,
                const bool followLink = true
            ) const = 0;

            //- Copy, recursively if necessary, the source to the destination
            virtual bool cp
            (
                const fileName& src,
                const fileName& dst,
                const bool followLink = true
            ) const = 0;

            //- Create a softlink. dst should not exist. Returns true if
            //  successful.
            virtual bool ln(const fileName& src, const fileName& dst) const = 0;

            //- Rename src to dst
            virtual bool mv
            (
                const fileName& src,
                const fileName& dst,
                const bool followLink = false
            ) const = 0;

            //- Rename to a corresponding backup file
            //  If the backup file already exists, attempt with
            //  "01" .. "99" suffix
            virtual bool mvBak
            (
                const fileName&,
                const std::string& ext = "bak"
            ) const = 0;

            //- Remove a file, returning true if successful otherwise false
            virtual bool rm(const fileName&) const = 0;

            //- Remove a directory and its contents
            //  \param silent do not report missing directory
            virtual bool rmDir
            (
                const fileName& dir,
                const bool silent = false
            ) const = 0;

//            //- Open a shared library. Return handle to library. Print error
//            //  message if library cannot be loaded (check = true)
//            virtual void* dlOpen
//            (
//                const fileName& lib,
//                const bool check = true
//            ) const = 0;


        // (reg)IOobject functionality

            //- Generate disk file name for object. Opposite of filePath.
            //  Optional wanted typeName.
            virtual fileName objectPath
            (
                const IOobject& io,
                const word& typeName
            ) const;

            //- Search for an object. checkGlobal : also check undecomposed case
            //  Optional wanted typeName.
            virtual fileName filePath
            (
                const bool checkGlobal,
                const IOobject&,
                const word& typeName,
                const bool search = true
            ) const = 0;

            //- Search for a directory. checkGlobal : also check undecomposed
            //  case
            virtual fileName dirPath
            (
                const bool checkGlobal,
                const IOobject& io,
                const bool search = true
            ) const = 0;

            //- Search directory for objects. Used in IOobjectList.
            virtual fileNameList readObjects
            (
                const objectRegistry& db,
                const fileName& instance,
                const fileName& local,
                word& newInstance
            ) const;

            //- Read object header from supplied file
            virtual bool readHeader
            (
                IOobject&,
                const fileName&,
                const word& typeName
            ) const = 0;

            //- Reads header for regIOobject and returns an ISstream
            //  to read the contents.
            virtual autoPtr<ISstream> readStream
            (
                regIOobject&,
                const fileName&,
                const word& typeName,
                const bool valid = true
            ) const = 0;

            //- Top-level read
            virtual bool read
            (
                regIOobject&,
                const bool masterOnly,
                const IOstreamOption::streamFormat format,
                const word& typeName
            ) const = 0;

            //- Writes a regIOobject (so header, contents and divider).
            //  Returns success state. Default action is to write to
            //  the objectPath using writeData. If !valid the
            //  file does not need to be written (this is used e.g. to
            //  suppress empty local lagrangian data)
            virtual bool writeObject
            (
                const regIOobject& io,
                IOstreamOption streamOpt = IOstreamOption(),
                const bool valid = true
            ) const;


        // Filename (not IOobject) operations

            //- Search for a file or directory. Use IOobject version in
            //  preference
            virtual fileName filePath(const fileName&) const;

            //- Generate an ISstream that reads a file
            virtual autoPtr<ISstream> NewIFstream(const fileName&) const = 0;

            //- Generate an OSstream that writes a file
            virtual autoPtr<OSstream> NewOFstream
            (
                const fileName& pathname,
                IOstreamOption streamOpt = IOstreamOption(),
                const bool valid = true
            ) const = 0;


        // File modification checking

            //- Add watching of a file. Returns handle
            virtual label addWatch(const fileName&) const;

            //- Remove watch on a file (using handle)
            virtual bool removeWatch(const label) const;

            //- Find index (or -1) of file in list of handles
            virtual label findWatch
            (
                const labelList& watchIndices,
                const fileName&
            ) const;

            //- Helper: add watches for list of regIOobjects
            virtual void addWatches(regIOobject&, const fileNameList&) const;

            //- Get name of file being watched (using handle)
            virtual fileName getFile(const label) const;

            //- Update state of all files
            virtual void updateStates
            (
                const bool masterOnly,
                const bool syncPar
            ) const;

            //- Get current state of file (using handle)
            virtual fileMonitor::fileState getState(const label) const;

            //- Set current state of file (using handle) to unmodified
            virtual void setUnmodified(const label) const;


        // Other

            //- Actual name of processors dir (for use in mode PROCOBJECT,
            //  PROCINSTANCE)
            virtual word processorsDir(const IOobject& io) const
            {
                return processorsBaseDir;
            }

            //- Actual name of processors dir (for use in mode PROCOBJECT,
            //  PROCINSTANCE)
            virtual word processorsDir(const fileName&) const
            {
                return processorsBaseDir;
            }

            //- Set number of processor directories/results. Only used in
            //  decomposePar
            virtual void setNProcs(const label nProcs);

            //- Get number of processor directories/results. Used for e.g.
            //  reconstructPar, argList checking
            virtual label nProcs
            (
                const fileName& dir,
                const fileName& local = ""
            ) const;

            //- Get sorted list of times
            virtual instantList findTimes(const fileName&, const word&) const;

            //- Find instance where IOobject is. Fails if cannot be found
            //  and readOpt() is MUST_READ/MUST_READ_IF_MODIFIED. Otherwise
            //  returns stopInstance.
            virtual IOobject findInstance
            (
                const IOobject& io,
                const scalar startValue,
                const word& stopInstance
            ) const;

            //- Callback for time change
            virtual void setTime(const Time&) const
            {}

            //- Forcibly wait until all output done. Flush any cached data
            virtual void flush() const;

            //- Generate path (like io.path) from root+casename with any
            //  'processorXXX' replaced by procDir (usually 'processsors')
            fileName processorsCasePath
            (
                const IOobject&,
                const word& procDir
            ) const;

            //- Generate path (like io.path) with provided instance and any
            //  'processorXXX' replaced by procDir (usually 'processsors')
            fileName processorsPath
            (
                const IOobject&,
                const word& instance,
                const word& procDir
            ) const;

            //- Operating on fileName: replace processorXXX with procDir
            fileName processorsPath(const fileName&, const word& procDir) const;

            //- Split objectPath into part before 'processor' and part after.
            //
            //  Returns -1 or processor number and optionally number
            //  of processors. Use with care.
            //  - path/"processor"+Foam::name(proci)/local reconstructs input
            //  - path/"processors"+Foam::name(nProcs)/local reconstructs
            //    collated processors equivalence
            static label splitProcessorPath
            (
                const fileName& objectPath,
                fileName& path,
                fileName& procDir,
                fileName& local,
                procRangeType& group,
                label& nProcs
            );

            //- Detect processor number from '/aa/bb/processorDDD/cc'
            static label detectProcessorPath(const fileName& objPath);
};


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// Note: defined in fileOperation.C

//- Get current file handler
const fileOperation& fileHandler();

//- Replace, reset file handler.
//  \return old handler on change, null otherwise
autoPtr<fileOperation> fileHandler(autoPtr<fileOperation>&& newHandler);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
