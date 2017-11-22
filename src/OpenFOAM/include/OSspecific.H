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

InNamespace
    Foam

Description
    Functions used by OpenFOAM that are specific to POSIX compliant
    operating systems and need to be replaced or emulated on other systems.

SourceFiles
    POSIX.C

\*---------------------------------------------------------------------------*/

#ifndef OSspecific_H
#define OSspecific_H

#include "fileNameList.H"
#include "stringList.H"

#include <sys/types.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class CStringList;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Return the PID of this process
pid_t pid();

//- Return the parent PID of this process
pid_t ppid();

//- Return the group PID of this process
pid_t pgid();

//- True if environment variable of given name is defined.
//  Using an empty name is a no-op and always returns false.
bool env(const std::string& envName);

//- Get environment value for given envName.
//  Return string() if the environment is undefined or envName is empty.
string getEnv(const std::string& envName);

//- Set an environment variable, return true on success.
//  Using an empty name is a no-op and always returns false.
bool setEnv(const word& name, const std::string& value, const bool overwrite);

//- Return the system's host name, as per hostname(1)
//  Optionally with the full name (as per the '-f' option)
string hostName(const bool full=false);

//- Return the system's domain name, as per hostname(1) with the '-d' option
string domainName();

//- Return the user's login name
string userName();

//- Is user administrator
bool isAdministrator();

//- Return home directory path name for the current user
fileName home();

//- Return home directory path name for a particular user
fileName home(const std::string& userName);

//- Return current working directory path name
fileName cwd();

//- Change current directory to the one specified and return true on success.
//  Using an empty name is a no-op and always returns false.
bool chDir(const fileName& dir);

//- Make a directory and return an error if it could not be created
//  and does not already exist.
//  Using an empty pathName is a no-op and always returns false.
bool mkDir(const fileName& pathName, mode_t mode=0777);

//- Set the file/directory mode, return true on success.
//  Using an empty name is a no-op and always returns false.
bool chMod(const fileName& name, const mode_t mode);

//- Return the file mode, normally following symbolic links
//  Using an empty name is a no-op and always returns 0.
mode_t mode(const fileName& name, const bool followLink=true);

//- Return the file type: DIRECTORY or FILE, normally following symbolic links
//  Using an empty name is a no-op and always returns UNDEFINED.
fileName::Type type(const fileName& name, const bool followLink=true);

//- Does the name exist (as DIRECTORY or FILE) in the file system?
//  Optionally enable/disable check for gzip file.
//  Using an empty name is a no-op and always returns false.
bool exists
(
    const fileName& name,
    const bool checkGzip=true,
    const bool followLink=true
);

//- Does the name exist as a DIRECTORY in the file system?
//  Using an empty name is a no-op and always returns false.
bool isDir(const fileName& name, const bool followLink=true);

//- Does the name exist as a FILE in the file system?
//  Optionally enable/disable check for gzip file.
//  Using an empty name is a no-op and always returns false.
bool isFile
(
    const fileName& name,
    const bool checkGzip=true,
    const bool followLink=true
);

//- Return size of file or -1 on failure (normally follows symbolic links).
//  Using an empty name is a no-op and always returns -1.
off_t fileSize(const fileName& name, const bool followLink=true);

//- Return time of last file modification (normally follows symbolic links).
//  Using an empty name is a no-op and always returns 0.
time_t lastModified(const fileName& name, const bool followLink=true);

//- Return time of last file modification
//  Using an empty name is a no-op and always returns 0.
double highResLastModified(const fileName&, const bool followLink = true);

//- Read a directory and return the entries as a fileName List.
//  Using an empty directory name returns an empty list.
fileNameList readDir
(
    const fileName& directory,
    const fileName::Type type=fileName::FILE,
    const bool filtergz=true,
    const bool followLink=true
);

//- Copy the source to the destination (recursively if necessary).
//  An empty source name is a no-op and always returns false.
bool cp(const fileName& src, const fileName& dst, const bool followLink=true);

//- Create a softlink. dst should not exist. Returns true if successful.
//  An empty source or destination name is a no-op that always returns false,
//  but also produces a warning.
bool ln(const fileName& src, const fileName& dst);

//- Rename src to dst.
//  An empty source or destination name is a no-op that always returns false.
bool mv
(
    const fileName& src,
    const fileName& dst,
    const bool followLink=false
);

//- Rename to a corresponding backup file
//  If the backup file already exists, attempt with "01" .. "99" suffix
//  An empty name or extension is a no-op that always returns false.
bool mvBak(const fileName& src, const std::string& ext = "bak");

//- Remove a file (or its gz equivalent), returning true if successful.
//  An empty name is a no-op that always returns false.
bool rm(const fileName& file);

//- Remove a dirctory and its contents (optionally silencing warnings)
//  An empty directory name is a no-op that always returns false,
//  but also produces a warning.
bool rmDir(const fileName& directory, const bool silent=false);

//- Sleep for the specified number of seconds
unsigned int sleep(const unsigned int sec);

//- Close file descriptor
void fdClose(const int fd);

//- Check if machine is up by pinging given port
bool ping(const std::string& destName, const label port, const label timeOut);

//- Check if machine is up by pinging port 22 (ssh) and 222 (rsh)
bool ping(const std::string& host, const label timeOut=10);

//- Execute the specified command via the shell.
//  Uses vfork/execl internally.
//  Where possible, use the list version instead.
//
//  \param background return immediately to parent process instead of waiting
//      for the child. Can be used (with moderation) to create background
//      processes.
//
//  \note treats an empty command as a successful no-op.
int system(const std::string& command, const bool background=false);

//- Execute the specified command with arguments.
//  Uses vfork/execvp internally
//
//  \param background return immediately to parent process instead of waiting
//      for the child. Can be used (with moderation) to create background
//      processes.
//
//  \note treats an empty command as a successful no-op.
int system(const UList<string>& command, const bool background=false);

//- Execute the specified command with arguments.
//  Uses vfork/execvp internally
//
//  \param background return immediately to parent process instead of waiting
//      for the child. Can be used (with moderation) to create background
//      processes.
//
//  \note treats an empty command as a successful no-op.
int system(const CStringList& command, const bool background=false);

//- Open a shared library and return handle to library.
//  Print error message if library cannot be loaded (suppress with check=true)
void* dlOpen(const fileName& lib, const bool check = true);

//- Close a dlopened library using handle. Return true if successful
bool dlClose(void* handle);

//- Lookup a symbol in a dlopened library using handle to library
void* dlSym(void* handle, const std::string& symbol);

//- Report if symbol in a dlopened library could be found.
//  Using a null handle or an empty symbol name is a no-op and always
//  returns false.
bool dlSymFound(void* handle, const std::string& symbol);

//- Return all loaded libraries
fileNameList dlLoaded();


// Low level random numbers. Use Random class instead.

//- Seed random number generator.
void osRandomSeed(const label seed);

//- Return random integer (uniform distribution between 0 and 2^31)
label osRandomInteger();

//- Return random double precision (uniform distribution between 0 and 1)
scalar osRandomDouble();

//- Return the size of the buffer for re-entrant random number generator
label osRandomBufferSize();

//- Seed random number generator.
void osRandomSeed(const label seed, List<char>& buffer);

//- Return random integer (uniform distribution between 0 and 2^31)
label osRandomInteger(List<char>& buffer);

//- Return random double precision (uniform distribution between 0 and 1)
scalar osRandomDouble(List<char>& buffer);


// Thread handling

//- Allocate a thread
label allocateThread();

//- Start a thread
void createThread(const label, void *(*start_routine) (void *), void *arg);

//- Wait for thread
void joinThread(const label);

//- Delete a thread
void freeThread(const label);

//- Allocate a mutex variable
label allocateMutex();

//- Lock a mutex variable
void lockMutex(const label);

//- Unlock a mutex variable
void unlockMutex(const label);

//- Free a mutex variable
void freeMutex(const label);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
