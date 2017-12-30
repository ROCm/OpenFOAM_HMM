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

Description
    POSIX versions of the functions declared in OSspecific.H

\*---------------------------------------------------------------------------*/

#ifdef solarisGcc
    #define _SYS_VNODE_H
#endif

#include "OSspecific.H"
#include "POSIX.H"
#include "foamVersion.H"
#include "fileName.H"
#include "fileStat.H"
#include "timer.H"
#include "IFstream.H"
#include "DynamicList.H"
#include "CStringList.H"
#include "SubList.H"
#include "IOstreams.H"
#include "Pstream.H"

#include <fstream>
#include <cstdlib>
#include <cctype>

#include <stdio.h>
#include <unistd.h>
#include <dirent.h>
#include <pwd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <sys/socket.h>
#include <netdb.h>
#include <dlfcn.h>
#include <link.h>

#include <netinet/in.h>
#ifdef USE_RANDOM
    #include <climits>
    #if INT_MAX    != 2147483647
        #error "INT_MAX    != 2147483647"
        #error "The random number generator may not work!"
    #endif
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(POSIX, 0);
}


// * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * * //

//
//! \cond fileScope
//
// Return true if filename appears to be a backup file
//
static inline bool isBackupName(const Foam::fileName& name)
{
    if (name.empty())
    {
        return false;
    }
    else if (name.back() == '~')
    {
        return true;
    }

    // Now check the extension
    const Foam::word ext = name.ext();
    if (ext.empty())
    {
        return false;
    }

    return
    (
        ext == "bak" || ext == "BAK"
     || ext == "old" || ext == "save"
    );
}


// Like fileName "/" global operator, but retain any invalid characters
static inline Foam::fileName fileNameConcat
(
    const std::string& a,
    const std::string& b
)
{
    if (a.size())
    {
        if (b.size())
        {
            // Two non-empty strings: can concatenate
            return Foam::fileName((a + '/' + b), false);
        }

        return Foam::fileName(a, false);
    }

    // Or, if the first string is empty

    if (b.size())
    {
        return Foam::fileName(b, false);
    }

    // Both strings are empty
    return Foam::fileName();
}
//! \endcond


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

pid_t Foam::pid()
{
    return ::getpid();
}


pid_t Foam::ppid()
{
    return ::getppid();
}


pid_t Foam::pgid()
{
    return ::getpgrp();
}


bool Foam::env(const std::string& envName)
{
    // An empty envName => always false
    return !envName.empty() && ::getenv(envName.c_str()) != nullptr;
}


Foam::string Foam::getEnv(const std::string& envName)
{
    // Ignore an empty envName => always ""
    char* env = envName.empty() ? nullptr : ::getenv(envName.c_str());

    if (env)
    {
        return string(env);
    }

    // Return null-constructed string rather than string::null
    // to avoid cyclic dependencies in the construction of globals
    return string();
}


bool Foam::setEnv
(
    const word& envName,
    const std::string& value,
    const bool overwrite
)
{
    // Ignore an empty envName => always false
    return
    (
        !envName.empty()
     && ::setenv(envName.c_str(), value.c_str(), overwrite) == 0
    );
}


Foam::string Foam::hostName(bool full)
{
    char buf[128];
    ::gethostname(buf, sizeof(buf));

    // implementation as per hostname from net-tools
    if (full)
    {
        struct hostent *hp = ::gethostbyname(buf);
        if (hp)
        {
            return hp->h_name;
        }
    }

    return buf;
}


Foam::string Foam::domainName()
{
    char buf[128];
    ::gethostname(buf, sizeof(buf));

    // implementation as per hostname from net-tools
    struct hostent *hp = ::gethostbyname(buf);
    if (hp)
    {
        char *p = ::strchr(hp->h_name, '.');
        if (p)
        {
            ++p;
            return p;
        }
    }

    return string::null;
}


Foam::string Foam::userName()
{
    struct passwd* pw = ::getpwuid(::getuid());
    if (pw != nullptr)
    {
        return pw->pw_name;
    }

    return string();
}


bool Foam::isAdministrator()
{
    return (::geteuid() == 0);
}


Foam::fileName Foam::home()
{
    char* env = ::getenv("HOME");
    if (env)
    {
        return fileName(env);
    }

    struct passwd* pw = ::getpwuid(::getuid());
    if (pw)
    {
        return pw->pw_dir;
    }

    return fileName();
}


Foam::fileName Foam::home(const std::string& userName)
{
    // An empty userName => same as home()
    if (userName.empty())
    {
        return Foam::home();
    }

    struct passwd* pw = ::getpwnam(userName.c_str());
    if (pw)
    {
        return pw->pw_dir;
    }

    return fileName();
}


Foam::fileName Foam::cwd()
{
    label pathLengthLimit = POSIX::pathLengthChunk;
    List<char> path(pathLengthLimit);

    // Resize path if getcwd fails with an ERANGE error
    while (pathLengthLimit == path.size())
    {
        if (::getcwd(path.data(), path.size()))
        {
            return path.data();
        }
        else if (errno == ERANGE)
        {
            // Increment path length up to the pathLengthMax limit
            if
            (
                (pathLengthLimit += POSIX::pathLengthChunk)
             >= POSIX::pathLengthMax
            )
            {
                FatalErrorInFunction
                    << "Attempt to increase path length beyond limit of "
                    << POSIX::pathLengthMax
                    << exit(FatalError);
            }

            path.setSize(pathLengthLimit);
        }
        else
        {
            break;
        }
    }

    FatalErrorInFunction
        << "Couldn't get the current working directory"
        << exit(FatalError);

    return fileName::null;
}


bool Foam::chDir(const fileName& dir)
{
    // Ignore an empty dir name => always false
    return !dir.empty() && ::chdir(dir.c_str()) == 0;
}


bool Foam::mkDir(const fileName& pathName, mode_t mode)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : pathName:" << pathName << " mode:" << mode
            << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // empty names are meaningless
    if (pathName.empty())
    {
        return false;
    }

    // Construct path directory if does not exist
    if (::mkdir(pathName.c_str(), mode) == 0)
    {
        // Directory made OK so return true
        return true;
    }
    else
    {
        switch (errno)
        {
            case EPERM:
            {
                FatalErrorInFunction
                    << "The filesystem containing " << pathName
                    << " does not support the creation of directories."
                    << exit(FatalError);

                return false;
            }

            case EEXIST:
            {
                // Directory already exists so simply return true
                return true;
            }

            case EFAULT:
            {
                FatalErrorInFunction
                    << "" << pathName
                    << " points outside your accessible address space."
                    << exit(FatalError);

                return false;
            }

            case EACCES:
            {
                FatalErrorInFunction
                    << "The parent directory does not allow write "
                       "permission to the process,"<< nl
                    << "or one of the directories in " << pathName
                    << " did not allow search (execute) permission."
                    << exit(FatalError);

                return false;
            }

            case ENAMETOOLONG:
            {
                FatalErrorInFunction
                    << "" << pathName << " is too long."
                    << exit(FatalError);

                return false;
            }

            case ENOENT:
            {
                // Part of the path does not exist so try to create it
                if (pathName.path().size() && mkDir(pathName.path(), mode))
                {
                    return mkDir(pathName, mode);
                }
                else
                {
                    FatalErrorInFunction
                        << "Couldn't create directory " << pathName
                        << exit(FatalError);

                    return false;
                }
            }

            case ENOTDIR:
            {
                FatalErrorInFunction
                    << "A component used as a directory in " << pathName
                    << " is not, in fact, a directory."
                    << exit(FatalError);

                return false;
            }

            case ENOMEM:
            {
                FatalErrorInFunction
                    << "Insufficient kernel memory was available to make "
                       "directory " << pathName << '.'
                    << exit(FatalError);

                return false;
            }

            case EROFS:
            {
                FatalErrorInFunction
                    << "" << pathName
                    << " refers to a file on a read-only filesystem."
                    << exit(FatalError);

                return false;
            }

            case ELOOP:
            {
                FatalErrorInFunction
                    << "Too many symbolic links were encountered in resolving "
                    << pathName << '.'
                    << exit(FatalError);

                return false;
            }

            case ENOSPC:
            {
                FatalErrorInFunction
                    << "The device containing " << pathName
                    << " has no room for the new directory or "
                    << "the user's disk quota is exhausted."
                    << exit(FatalError);

                return false;
            }

            default:
            {
                FatalErrorInFunction
                    << "Couldn't create directory " << pathName
                    << exit(FatalError);

                return false;
            }
        }
    }
}


bool Foam::chMod(const fileName& name, const mode_t m)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Ignore an empty name => always false
    return !name.empty() && ::chmod(name.c_str(), m) == 0;
}


mode_t Foam::mode(const fileName& name, const bool followLink)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << endl;
    }

    // Ignore an empty name => always 0
    if (!name.empty())
    {
        fileStat fileStatus(name, followLink);
        if (fileStatus.isValid())
        {
            return fileStatus.status().st_mode;
        }
    }

    return 0;
}


Foam::fileName::Type Foam::type(const fileName& name, const bool followLink)
{
    // Ignore an empty name => always UNDEFINED
    if (name.empty())
    {
        return fileName::UNDEFINED;
    }

    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    mode_t m = mode(name, followLink);

    if (S_ISREG(m))
    {
        return fileName::FILE;
    }
    else if (S_ISLNK(m))
    {
        return fileName::LINK;
    }
    else if (S_ISDIR(m))
    {
        return fileName::DIRECTORY;
    }

    return fileName::UNDEFINED;
}


bool Foam::exists
(
    const fileName& name,
    const bool checkGzip,
    const bool followLink
)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << " checkGzip:" << checkGzip
            << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Ignore an empty name => always false
    return
    (
        !name.empty()
     && (mode(name, followLink) || isFile(name, checkGzip, followLink))
    );
}


bool Foam::isDir(const fileName& name, const bool followLink)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Ignore an empty name => always false
    return !name.empty() && S_ISDIR(mode(name, followLink));
}


bool Foam::isFile
(
    const fileName& name,
    const bool checkGzip,
    const bool followLink
)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << " checkGzip:" << checkGzip
            << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Ignore an empty name => always false
    return
    (
        !name.empty()
     && (
            S_ISREG(mode(name, followLink))
         || (checkGzip && S_ISREG(mode(name + ".gz", followLink)))
        )
    );
}


off_t Foam::fileSize(const fileName& name, const bool followLink)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Ignore an empty name
    if (!name.empty())
    {
        fileStat fileStatus(name, followLink);
        if (fileStatus.isValid())
        {
            return fileStatus.status().st_size;
        }
    }

    return -1;
}


time_t Foam::lastModified(const fileName& name, const bool followLink)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Ignore an empty name
    if (!name.empty())
    {
        fileStat fileStatus(name, followLink);
        if (fileStatus.isValid())
        {
            return fileStatus.status().st_mtime;
        }
    }

    return 0;
}


double Foam::highResLastModified(const fileName& name, const bool followLink)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Ignore an empty name
    if (!name.empty())
    {
        fileStat fileStatus(name);
        if (fileStatus.isValid())
        {
            return
                fileStatus.status().st_mtime
              + 1e-9*fileStatus.status().st_atim.tv_nsec;
        }
    }

    return 0;
}


Foam::fileNameList Foam::readDir
(
    const fileName& directory,
    const fileName::Type type,
    const bool filtergz,
    const bool followLink
)
{
    // Initial filename list size and the increment when resizing the list
    static const int maxNnames = 100;

    // Basic sanity: cannot strip '.gz' from directory names
    const bool stripgz = filtergz && (type != fileName::DIRECTORY);
    const word extgz("gz");

    fileNameList dirEntries;

    // Open directory and set the structure pointer
    // Do not attempt to open an empty directory name
    DIR *source;
    if
    (
        directory.empty()
     || (source = ::opendir(directory.c_str())) == nullptr
    )
    {
        if (POSIX::debug)
        {
            InfoInFunction
                << "cannot open directory " << directory << endl;
        }

        return dirEntries;
    }

    if (POSIX::debug)
    {
        // InfoInFunction
        Pout<< FUNCTION_NAME << " : reading directory " << directory << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    label nFailed = 0;     // Entries with invalid characters
    label nEntries = 0;    // Number of selected entries
    dirEntries.setSize(maxNnames);

    // Read and parse all the entries in the directory
    for (struct dirent *list; (list = ::readdir(source)) != nullptr; /*nil*/)
    {
        const std::string item(list->d_name);

        // Ignore files/directories beginning with "."
        // These are the ".", ".." directories and any hidden files/dirs
        if (item.empty() || item[0] == '.')
        {
            continue;
        }

        // Validate filename without spaces, quotes, etc in the name.
        // No duplicate slashes to strip - dirent will not have them anyhow.

        const fileName name(fileName::validate(item));
        if (name != item)
        {
            ++nFailed;
        }
        else if
        (
            (type == fileName::DIRECTORY)
         || (type == fileName::FILE && !isBackupName(name))
        )
        {
            if ((directory/name).type(followLink) == type)
            {
                if (nEntries >= dirEntries.size())
                {
                    dirEntries.setSize(dirEntries.size() + maxNnames);
                }

                if (stripgz && name.hasExt(extgz))
                {
                    dirEntries[nEntries++] = name.lessExt();
                }
                else
                {
                    dirEntries[nEntries++] = name;
                }
            }
        }
    }
    ::closedir(source);

    // Finalize the length of the entries list
    dirEntries.setSize(nEntries);

    if (nFailed && POSIX::debug)
    {
        std::cerr
            << "Foam::readDir() : reading directory " << directory << nl
            << nFailed << " entries with invalid characters in their name"
            << std::endl;
    }

    return dirEntries;
}


bool Foam::cp(const fileName& src, const fileName& dest, const bool followLink)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : src:" << src << " dest:" << dest << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Make sure source exists - this also handles an empty source name
    if (!exists(src))
    {
        return false;
    }

    const fileName::Type srcType = src.type(followLink);

    fileName destFile(dest);

    // Check type of source file.
    if (srcType == fileName::FILE)
    {
        // If dest is a directory, create the destination file name.
        if (destFile.type() == fileName::DIRECTORY)
        {
            destFile = destFile/src.name();
        }

        // Make sure the destination directory exists.
        if (!isDir(destFile.path()) && !mkDir(destFile.path()))
        {
            return false;
        }

        // Open and check streams.
        std::ifstream srcStream(src);
        if (!srcStream)
        {
            return false;
        }

        std::ofstream destStream(destFile);
        if (!destStream)
        {
            return false;
        }

        // Copy character data.
        char ch;
        while (srcStream.get(ch))
        {
            destStream.put(ch);
        }

        // Final check.
        if (!srcStream.eof() || !destStream)
        {
            return false;
        }
    }
    else if (srcType == fileName::LINK)
    {
        // If dest is a directory, create the destination file name.
        if (destFile.type() == fileName::DIRECTORY)
        {
            destFile = destFile/src.name();
        }

        // Make sure the destination directory exists.
        if (!isDir(destFile.path()) && !mkDir(destFile.path()))
        {
            return false;
        }

        ln(src, destFile);
    }
    else if (srcType == fileName::DIRECTORY)
    {
        // If dest is a directory, create the destination file name.
        if (destFile.type() == fileName::DIRECTORY)
        {
            destFile = destFile/src.components().last();
        }

        // Make sure the destination directory exists.
        if (!isDir(destFile) && !mkDir(destFile))
        {
            return false;
        }

        char* realSrcPath = realpath(src.c_str(), nullptr);
        char* realDestPath = realpath(destFile.c_str(), nullptr);
        const bool samePath = strcmp(realSrcPath, realDestPath) == 0;

        if (POSIX::debug && samePath)
        {
            InfoInFunction
                << "Attempt to copy " << realSrcPath << " to itself" << endl;
        }

        if (realSrcPath)
        {
            free(realSrcPath);
        }

        if (realDestPath)
        {
            free(realDestPath);
        }

        // Do not copy over self when src is actually a link to dest
        if (samePath)
        {
            return false;
        }

        // Copy files
        fileNameList files = readDir(src, fileName::FILE, false, followLink);
        for (const fileName& item : files)
        {
            if (POSIX::debug)
            {
                InfoInFunction
                    << "Copying : " << src/item
                    << " to " << destFile/item << endl;
            }

            // File to file.
            cp(src/item, destFile/item, followLink);
        }

        // Copy sub directories.
        fileNameList dirs = readDir
        (
            src,
            fileName::DIRECTORY,
            false,
            followLink
        );

        for (const fileName& item : dirs)
        {
            if (POSIX::debug)
            {
                InfoInFunction
                    << "Copying : " << src/item
                    << " to " << destFile << endl;
            }

            // Dir to Dir.
            cp(src/item, destFile, followLink);
        }
    }
    else
    {
        return false;
    }

    return true;
}


bool Foam::ln(const fileName& src, const fileName& dst)
{
    if (POSIX::debug)
    {
        //InfoInFunction
        Pout<< FUNCTION_NAME
            << " : Create softlink from : " << src << " to " << dst << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    if (src.empty())
    {
        WarningInFunction
            << "source name is empty: not linking." << endl;
        return false;
    }

    if (dst.empty())
    {
        WarningInFunction
            << "destination name is empty: not linking." << endl;
        return false;
    }

    if (exists(dst))
    {
        WarningInFunction
            << "destination " << dst << " already exists. Not linking."
            << endl;
        return false;
    }

    if (src.isAbsolute() && !exists(src))
    {
        WarningInFunction
            << "source " << src << " does not exist." << endl;
        return false;
    }

    if (::symlink(src.c_str(), dst.c_str()) == 0)
    {
        return true;
    }

    WarningInFunction
        << "symlink from " << src << " to " << dst << " failed." << endl;
    return false;
}


bool Foam::mv(const fileName& src, const fileName& dst, const bool followLink)
{
    if (POSIX::debug)
    {
        //InfoInFunction
        Pout<< FUNCTION_NAME << " : Move : " << src << " to " << dst << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Ignore an empty names => always false
    if (src.empty() || dst.empty())
    {
        return false;
    }

    if
    (
        dst.type() == fileName::DIRECTORY
     && src.type(followLink) != fileName::DIRECTORY
    )
    {
        const fileName dstName(dst/src.name());

        return ::rename(src.c_str(), dstName.c_str()) == 0;
    }
    else
    {
        return ::rename(src.c_str(), dst.c_str()) == 0;
    }
}


bool Foam::mvBak(const fileName& src, const std::string& ext)
{
    if (POSIX::debug)
    {
        //InfoInFunction
        Pout<< FUNCTION_NAME
            << " : moving : " << src << " to extension " << ext << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Ignore an empty name or extension => always false
    if (src.empty() || ext.empty())
    {
        return false;
    }

    if (exists(src, false))
    {
        const int maxIndex = 99;
        char index[3];

        for (int n = 0; n <= maxIndex; ++n)
        {
            fileName dstName(src + "." + ext);
            if (n)
            {
                ::sprintf(index, "%02d", n);
                dstName += index;
            }

            // avoid overwriting existing files, except for the last
            // possible index where we have no choice
            if (!exists(dstName, false) || n == maxIndex)
            {
                return ::rename(src.c_str(), dstName.c_str()) == 0;
            }
        }
    }

    // fallthrough: nothing to do
    return false;
}


bool Foam::rm(const fileName& file)
{
    if (POSIX::debug)
    {
        //InfoInFunction
        Pout<< FUNCTION_NAME << " : Removing : " << file << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Ignore an empty name => always false
    if (file.empty())
    {
        return false;
    }

    // Try returning plain file name; if not there, try with .gz
    if (::remove(file.c_str()) == 0)
    {
        return true;
    }
    else
    {
        return ::remove(string(file + ".gz").c_str()) == 0;
    }
}


bool Foam::rmDir(const fileName& directory, const bool silent)
{
    // Open directory and set the structure pointer
    // Do not attempt to open an empty directory name
    DIR *source;
    if
    (
        directory.empty()
     || (source = ::opendir(directory.c_str())) == nullptr
    )
    {
        if (!silent)
        {
            WarningInFunction
                << "cannot open directory " << directory << endl;
        }

        return false;
    }

    if (POSIX::debug)
    {
        //InfoInFunction
        Pout<< FUNCTION_NAME << " : removing directory " << directory << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Process each directory entry, counting any errors encountered
    label nErrors = 0;
    for (struct dirent *list; (list = ::readdir(source)) != nullptr; /*nil*/)
    {
        const std::string item(list->d_name);

        // Ignore "." and ".." directories
        if (item.empty() || item == "." || item == "..")
        {
            continue;
        }

        // Allow invalid characters (spaces, quotes, etc),
        // otherwise we cannot subdirs with these types of names.
        // -> const fileName path = directory/name; <-

        const fileName path(fileNameConcat(directory, item));

        if (path.type(false) == fileName::DIRECTORY)
        {
            if (!rmDir(path, true))  // Only report errors at the top-level
            {
                ++nErrors;
            }
        }
        else
        {
            if (!rm(path))
            {
                ++nErrors;
            }
        }
    }

    if (nErrors)
    {
        if (!silent)
        {
            WarningInFunction
                << "failed to remove directory " << directory << nl
                << "could not remove " << nErrors << " sub-entries" << endl;
        }
    }
    else
    {
        if (!rm(directory))
        {
            ++nErrors;
            if (!silent)
            {
                WarningInFunction
                    << "failed to remove directory " << directory << endl;
            }
        }
    }

    // clean up
    ::closedir(source);
    return !nErrors;
}


unsigned int Foam::sleep(const unsigned int sec)
{
    return ::sleep(sec);
}


void Foam::fdClose(const int fd)
{
    if (close(fd) != 0)
    {
        FatalErrorInFunction
            << "close error on " << fd << endl
            << abort(FatalError);
    }
}


bool Foam::ping
(
    const std::string& destName,
    const label destPort,
    const label timeOut
)
{
    struct hostent *hostPtr;
    volatile int sockfd;
    struct sockaddr_in destAddr;      // will hold the destination addr
    u_int addr;

    if ((hostPtr = ::gethostbyname(destName.c_str())) == nullptr)
    {
        FatalErrorInFunction
            << "gethostbyname error " << h_errno << " for host " << destName
            << abort(FatalError);
    }

    // Get first of the SLL of addresses
    addr = (reinterpret_cast<struct in_addr*>(*(hostPtr->h_addr_list)))->s_addr;

    // Allocate socket
    sockfd = ::socket(AF_INET, SOCK_STREAM, 0);
    if (sockfd < 0)
    {
        FatalErrorInFunction
            << "socket error"
            << abort(FatalError);
    }

    // Fill sockaddr_in structure with dest address and port
    memset(reinterpret_cast<char *>(&destAddr), '\0', sizeof(destAddr));
    destAddr.sin_family = AF_INET;
    destAddr.sin_port = htons(ushort(destPort));
    destAddr.sin_addr.s_addr = addr;


    timer myTimer(timeOut);

    if (timedOut(myTimer))
    {
        // Setjmp from timer jumps back to here
        fdClose(sockfd);
        return false;
    }

    if
    (
        ::connect
        (
            sockfd,
            reinterpret_cast<struct sockaddr*>(&destAddr),
            sizeof(struct sockaddr)
        ) != 0
    )
    {
        // Connection refused. Check if network was actually used or not.

        int connectErr = errno;

        fdClose(sockfd);

        if (connectErr == ECONNREFUSED)
        {
            return true;
        }
        //perror("connect");

        return false;
    }

    fdClose(sockfd);

    return true;
}


bool Foam::ping(const std::string& host, const label timeOut)
{
    return ping(host, 222, timeOut) || ping(host, 22, timeOut);
}


namespace Foam
{
//! \cond fileScope
static int waitpid(const pid_t pid)
{
    // child status, return code from the exec etc.
    int status = 0;

    // in parent - blocking wait
    // modest treatment of signals (in child)
    // treat 'stopped' like exit (suspend/continue)
    while (true)
    {
        pid_t wpid = ::waitpid(pid, &status, WUNTRACED);

        if (wpid == -1)
        {
            FatalErrorInFunction
                << "some error occurred in child"
                << exit(FatalError);
            break;
        }

        if (WIFEXITED(status))
        {
            // child exited, get its return status
            return WEXITSTATUS(status);
        }

        if (WIFSIGNALED(status))
        {
            // child terminated by some signal
            return WTERMSIG(status);
        }

        if (WIFSTOPPED(status))
        {
            // child stopped by some signal
            return WSTOPSIG(status);
        }

        FatalErrorInFunction
            << "programming error, status from waitpid() not handled: "
            << status
            << exit(FatalError);
    }

    return -1;  // should not happen
}
//! \endcond
}


int Foam::system
(
    const std::string& command,
    const bool background
)
{
    if (command.empty())
    {
        // Treat an empty command as a successful no-op.
        // From 'man sh' POSIX (man sh):
        //   "If the command_string operand is an empty string,
        //    sh shall exit with a zero exit status."
        return 0;
    }

    pid_t child_pid = ::vfork();   // NB: vfork, not fork!
    if (child_pid == -1)
    {
        FatalErrorInFunction
            << "vfork() failed for system command " << command
            << exit(FatalError);
    }

    if (child_pid == 0)
    {
        // in child

        // execl uses the current environ
        (void) ::execl
        (
            "/bin/sh",          // Path of the shell
            "sh",               // Command-name (name for the shell)
            "-c",               // Read commands from command_string operand
            command.c_str(),    // Command string
            reinterpret_cast<char*>(0)
        );

        // obviously failed, since exec should not return at all
        FatalErrorInFunction
            << "exec failed: " << command
            << exit(FatalError);
    }


    // In parent:

    if (background)
    {
        // Started as background process
        return 0;
    }

    // blocking wait for the child
    return waitpid(child_pid);
}


int Foam::system
(
    const CStringList& command,
    const bool background
)
{
    const int argc = command.size();

    if (!argc)
    {
        // Treat an empty command as a successful no-op.
        // For consistency with POSIX (man sh) behaviour for (sh -c command),
        // which is what is mostly being replicated here.
        return 0;
    }

    // NB: use vfork, not fork!
    // vfork behaves more like a thread and avoids copy-on-write problems
    // triggered by fork.
    // The normal system() command has a fork buried in it that causes
    // issues with infiniband and openmpi etc.
    pid_t child_pid = ::vfork();
    if (child_pid == -1)
    {
        FatalErrorInFunction
            << "vfork() failed for system command " << command[0]
            << exit(FatalError);
    }

    if (child_pid == 0)
    {
        // In child:
        // Need command and arguments separately.
        // args is a nullptr-terminated list of c-strings

        // execvp uses the current environ
        (void) ::execvp(command[0], command.strings(1));

        // obviously failed, since exec should not return at all
        FatalErrorInFunction
            << "exec(" << command[0] << ", ...) failed"
            << exit(FatalError);
    }


    // In parent:

    if (background)
    {
        // Started as background process
        return 0;
    }

    // blocking wait for the child
    return waitpid(child_pid);
}


int Foam::system
(
    const Foam::UList<Foam::string>& command,
    const bool background
)
{
    // In the future simply call the CStringList version:
    //
    //     const CStringList cmd(command);
    //     return Foam::system(cmd, background);

    const int argc = command.size();

    if (!argc)
    {
        // Treat an empty command as a successful no-op.
        // For consistency with POSIX (man sh) behaviour for (sh -c command),
        // which is what is mostly being replicated here.
        return 0;
    }

    // NB: use vfork, not fork!
    // vfork behaves more like a thread and avoids copy-on-write problems
    // triggered by fork.
    // The normal system() command has a fork buried in it that causes
    // issues with infiniband and openmpi etc.
    pid_t child_pid = ::vfork();
    if (child_pid == -1)
    {
        FatalErrorInFunction
            << "vfork() failed for system command " << command[0]
            << exit(FatalError);
    }

    if (child_pid == 0)
    {
        // In child:
        // Need command and arguments separately.
        // args is a nullptr-terminated list of c-strings

        CStringList args(SubList<string>(command, 0));
        if (argc > 1)
        {
            args.reset(SubList<string>(command, argc-1, 1));
        }

        // execvp uses the current environ
        (void) ::execvp(command[0].c_str(), args.strings());

        // obviously failed, since exec should not return at all
        FatalErrorInFunction
            << "exec(" << command[0] << ", ...) failed"
            << exit(FatalError);
    }


    // In parent:

    if (background)
    {
        // Started as background process
        return 0;
    }

    // blocking wait for the child
    return waitpid(child_pid);
}


void* Foam::dlOpen(const fileName& lib, const bool check)
{
    if (POSIX::debug)
    {
        std::cout<< "dlOpen(const fileName&)"
            << " : dlopen of " << lib << std::endl;
    }
    void* handle = ::dlopen(lib.c_str(), RTLD_LAZY|RTLD_GLOBAL);

    if (!handle && check)
    {
        WarningInFunction
            << "dlopen error : " << ::dlerror()
            << endl;
    }

    if (POSIX::debug)
    {
        std::cout
            << "dlOpen(const fileName&)"
            << " : dlopen of " << lib
            << " handle " << handle << std::endl;
    }

    return handle;
}


bool Foam::dlClose(void* handle)
{
    if (POSIX::debug)
    {
        std::cout
            << "dlClose(void*)"
            << " : dlclose of handle " << handle << std::endl;
    }
    return ::dlclose(handle) == 0;
}


void* Foam::dlSym(void* handle, const std::string& symbol)
{
    if (POSIX::debug)
    {
        std::cout
            << "dlSym(void*, const std::string&)"
            << " : dlsym of " << symbol << std::endl;
    }
    // clear any old errors - see manpage dlopen
    (void) ::dlerror();

    // get address of symbol
    void* fun = ::dlsym(handle, symbol.c_str());

    // find error (if any)
    char *error = ::dlerror();

    if (error)
    {
        WarningInFunction
            << "Cannot lookup symbol " << symbol << " : " << error
            << endl;
    }

    return fun;
}


bool Foam::dlSymFound(void* handle, const std::string& symbol)
{
    if (handle && !symbol.empty())
    {
        if (POSIX::debug)
        {
            std::cout
                << "dlSymFound(void*, const std::string&)"
                << " : dlsym of " << symbol << std::endl;
        }

        // clear any old errors - see manpage dlopen
        (void) ::dlerror();

        // get address of symbol
        (void) ::dlsym(handle, symbol.c_str());

        // symbol can be found if there was no error
        return !::dlerror();
    }

    return false;
}


static int collectLibsCallback
(
    struct dl_phdr_info *info,
    size_t size,
    void *data
)
{
    Foam::DynamicList<Foam::fileName>* ptr =
        reinterpret_cast<Foam::DynamicList<Foam::fileName>*>(data);
    ptr->append(info->dlpi_name);
    return 0;
}


Foam::fileNameList Foam::dlLoaded()
{
    DynamicList<fileName> libs;
    dl_iterate_phdr(collectLibsCallback, &libs);
    if (POSIX::debug)
    {
        std::cout
            << "dlLoaded()"
            << " : determined loaded libraries :" << libs.size() << std::endl;
    }
    return libs;
}

Foam::label Foam::osRandomBufferSize()
{
    #ifdef USE_RANDOM
    return sizeof(random_data);
    #else
    return sizeof(drand48_data);
    #endif
}


void Foam::osRandomSeed(const label seed)
{
    #ifdef USE_RANDOM
    srandom((unsigned int)seed);
    #else
    srand48(seed);
    #endif
}


Foam::label Foam::osRandomInteger()
{
    #ifdef USE_RANDOM
    return random();
    #else
    return lrand48();
    #endif
}


Foam::scalar Foam::osRandomDouble()
{
    #ifdef USE_RANDOM
    return (scalar)random()/INT_MAX;
    #else
    return drand48();
    #endif
}


void Foam::osRandomSeed(const label seed, List<char>& buffer)
{
    #ifdef USE_RANDOM
    srandom_r((unsigned int)seed, reinterpret_cast<random_data*>(buffer.begin()));
    #else
    srand48_r(seed, reinterpret_cast<drand48_data*>(buffer.begin()));
    #endif
}


Foam::label Foam::osRandomInteger(List<char>& buffer)
{
    #ifdef USE_RANDOM
    int32_t result;
    random_r(reinterpret_cast<random_data*>(buffer.begin()), &result);
    return result;
    #else
    long result;
    lrand48_r(reinterpret_cast<drand48_data*>(buffer.begin()), &result);
    return result;
    #endif
}


Foam::scalar Foam::osRandomDouble(List<char>& buffer)
{
    #ifdef USE_RANDOM
    int32_t result;
    random_r(reinterpret_cast<random_data*>(buffer.begin()), &result);
    return (scalar)result/INT_MAX;
    #else
    double result;
    drand48_r(reinterpret_cast<drand48_data*>(buffer.begin()), &result);
    return result;
    #endif

}


static Foam::DynamicList<Foam::autoPtr<pthread_t>> threads_;
static Foam::DynamicList<Foam::autoPtr<pthread_mutex_t>> mutexes_;

Foam::label Foam::allocateThread()
{
    forAll(threads_, i)
    {
        if (!threads_[i].valid())
        {
            if (POSIX::debug)
            {
                Pout<< "allocateThread : reusing index:" << i << endl;
            }
            // Reuse entry
            threads_[i].reset(new pthread_t());
            return i;
        }
    }

    const label index = threads_.size();
    if (POSIX::debug)
    {
        Pout<< "allocateThread : new index:" << index << endl;
    }
    threads_.append(autoPtr<pthread_t>(new pthread_t()));

    return index;
}


void Foam::createThread
(
    const label index,
    void *(*start_routine) (void *),
    void *arg
)
{
    if (POSIX::debug)
    {
        Pout<< "createThread : index:" << index << endl;
    }
    if (pthread_create(&threads_[index](), nullptr, start_routine, arg))
    {
        FatalErrorInFunction
            << "Failed starting thread " << index << exit(FatalError);
    }
}


void Foam::joinThread(const label index)
{
    if (POSIX::debug)
    {
        Pout<< "freeThread : join:" << index << endl;
    }
    if (pthread_join(threads_[index](), nullptr))
    {
        FatalErrorInFunction << "Failed freeing thread " << index
            << exit(FatalError);
    }
}


void Foam::freeThread(const label index)
{
    if (POSIX::debug)
    {
        Pout<< "freeThread : index:" << index << endl;
    }
    threads_[index].clear();
}


Foam::label Foam::allocateMutex()
{
    forAll(mutexes_, i)
    {
        if (!mutexes_[i].valid())
        {
            if (POSIX::debug)
            {
                Pout<< "allocateMutex : reusing index:" << i << endl;
            }
            // Reuse entry
            mutexes_[i].reset(new pthread_mutex_t());
            return i;
        }
    }

    const label index = mutexes_.size();

    if (POSIX::debug)
    {
        Pout<< "allocateMutex : new index:" << index << endl;
    }
    mutexes_.append(autoPtr<pthread_mutex_t>(new pthread_mutex_t()));
    return index;
}


void Foam::lockMutex(const label index)
{
    if (POSIX::debug)
    {
        Pout<< "lockMutex : index:" << index << endl;
    }
    if (pthread_mutex_lock(&mutexes_[index]()))
    {
        FatalErrorInFunction << "Failed locking mutex " << index
            << exit(FatalError);
    }
}


void Foam::unlockMutex(const label index)
{
    if (POSIX::debug)
    {
        Pout<< "unlockMutex : index:" << index << endl;
    }
    if (pthread_mutex_unlock(&mutexes_[index]()))
    {
        FatalErrorInFunction << "Failed unlocking mutex " << index
            << exit(FatalError);
    }
}


void Foam::freeMutex(const label index)
{
    if (POSIX::debug)
    {
        Pout<< "freeMutex : index:" << index << endl;
    }
    mutexes_[index].clear();
}


// ************************************************************************* //
