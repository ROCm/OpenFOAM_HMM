/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "fstreamPointer.H"
#include "OCountStream.H"
#include "OSspecific.H"

// HAVE_LIBZ defined externally
// #define HAVE_LIBZ

#ifdef HAVE_LIBZ
#include "gzstream.h"
#endif /* HAVE_LIBZ */

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    static inline void removeConflictingFiles
    (
        const fileName& otherName,
        const bool append,
        const fileName& targetName
    )
    {
        // Remove other (compressed/uncompressed) version

        const fileName::Type pathType = Foam::type(otherName, false);

        if (pathType == fileName::FILE || pathType == fileName::LINK)
        {
            Foam::rm(otherName);
        }

        // Disallow writing into symlinked files.
        // Eg, avoid problems with symlinked initial fields

        if (!append && Foam::type(targetName, false) == fileName::LINK)
        {
            Foam::rm(targetName);
        }
    }
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::ifstreamPointer::supports_gz()
{
    #ifdef HAVE_LIBZ
    return true;
    #else
    return false;
    #endif
}


bool Foam::ofstreamPointer::supports_gz()
{
    #ifdef HAVE_LIBZ
    return true;
    #else
    return false;
    #endif
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ifstreamPointer::ifstreamPointer
(
    const fileName& pathname
)
:
    ptr_(nullptr)
{
    const std::ios_base::openmode mode
    (
        std::ios_base::in | std::ios_base::binary
    );

    ptr_.reset(new std::ifstream(pathname, mode));

    if (!ptr_->good())
    {
        // Try compressed version instead

        const fileName pathname_gz(pathname + ".gz");

        if (isFile(pathname_gz, false))
        {
            #ifdef HAVE_LIBZ

            ptr_.reset(new igzstream(pathname_gz, mode));

            #else /* HAVE_LIBZ */

            FatalError
                << "No read support for gz compressed files (libz)"
                << " : could use 'gunzip' from the command-line" << nl
                << "file: " << pathname_gz << endl
                << exit(FatalError);

            #endif /* HAVE_LIBZ */
        }
    }
}


Foam::ofstreamPointer::ofstreamPointer(std::nullptr_t)
:
    ptr_(new Foam::ocountstream)
{}


Foam::ofstreamPointer::ofstreamPointer
(
    const fileName& pathname,
    IOstreamOption::compressionType comp,
    const bool append
)
:
    ptr_(nullptr)
{
    std::ios_base::openmode mode
    (
        std::ios_base::out | std::ios_base::binary
    );

    if (append)
    {
        mode |= std::ios_base::app;
    }

    const fileName pathname_gz(pathname + ".gz");

    if (IOstreamOption::COMPRESSED == comp)
    {
        // Output compression requested

        #ifdef HAVE_LIBZ

        removeConflictingFiles(pathname, append, pathname_gz);
        ptr_.reset(new ogzstream(pathname_gz, mode));

        #else /* HAVE_LIBZ */

        comp = IOstreamOption::UNCOMPRESSED;

        Warning
            << nl
            << "No write support for gz compressed files (libz)"
            << " : downgraded to UNCOMPRESSED" << nl
            << "file: " << pathname_gz << endl;

        #endif /* HAVE_LIBZ */
    }

    if (IOstreamOption::UNCOMPRESSED == comp)
    {
        removeConflictingFiles(pathname_gz, append, pathname);
        ptr_.reset(new std::ofstream(pathname, mode));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ifstreamPointer::reopen_gz(const std::string& pathname_gz)
{
    #ifdef HAVE_LIBZ
    igzstream* gz = dynamic_cast<igzstream*>(ptr_.get());

    if (gz)
    {
        // Special treatment for gzstream
        gz->close();
        gz->clear();
        gz->open(pathname_gz);
    }
    #endif /* HAVE_LIBZ */
}


void Foam::ofstreamPointer::reopen_gz(const std::string& pathname_gz)
{
    #ifdef HAVE_LIBZ
    ogzstream* gz = dynamic_cast<ogzstream*>(ptr_.get());

    if (gz)
    {
        // Special treatment for gzstream
        gz->close();
        gz->clear();
        gz->open(pathname_gz);
    }
    #endif /* HAVE_LIBZ */
}


void Foam::ofstreamPointer::reopen(const std::string& pathname)
{
    std::ofstream* file = dynamic_cast<std::ofstream*>(ptr_.get());

    if (file)
    {
        if (file->is_open())
        {
            file->close();
        }
        file->clear();
        file->open(pathname);
    }
}


Foam::IOstreamOption::compressionType
Foam::ifstreamPointer::whichCompression() const
{
    #ifdef HAVE_LIBZ
    if (dynamic_cast<const igzstream*>(ptr_.get()))
    {
        return IOstreamOption::compressionType::COMPRESSED;
    }
    #endif /* HAVE_LIBZ */

    return IOstreamOption::compressionType::UNCOMPRESSED;
}


Foam::IOstreamOption::compressionType
Foam::ofstreamPointer::whichCompression() const
{
    #ifdef HAVE_LIBZ
    if (dynamic_cast<const ogzstream*>(ptr_.get()))
    {
        return IOstreamOption::compressionType::COMPRESSED;
    }
    #endif /* HAVE_LIBZ */

    return IOstreamOption::compressionType::UNCOMPRESSED;
}


// ************************************************************************* //
