/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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
    Foam::fileStat

Description
    Wrapper for stat() and lstat() system calls.

Warning
    On Linux (an maybe on others) a stat() of an nfs mounted (remote)
    file does never timeout and cannot be interrupted!
    So e.g. Foam::ping first and hope nfs is running.

SourceFiles
    fileStat.C

\*---------------------------------------------------------------------------*/

#ifndef fileStat_H
#define fileStat_H

#include <sys/stat.h>
#include <sys/types.h>

#include "label.H"
#include "fileName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class fileStat;

Istream& operator>>(Istream& is, fileStat& fs);
Ostream& operator<<(Ostream& os, const fileStat& fs);


/*---------------------------------------------------------------------------*\
                          Class fileStat Declaration
\*---------------------------------------------------------------------------*/

class fileStat
{
    // Private data

        struct stat status_;

        bool valid_;


public:

    // Constructors

        //- Empty constructor
        fileStat();

        //- Construct from components.
        //
        //  \param fName The file name or directory name to stat.
        //  \param followLink If it is a link, get the status of the source
        //      file/directory.
        //  \param maxTime The timeout value.
        //
        //  \note An empty filename is a no-op.
        fileStat
        (
            const char* fName,
            const bool followLink = true,
            const unsigned int maxTime = 0
        );

        //- Construct from components.
        //
        //  \param fName The file name or directory name to stat.
        //  \param followLink If it is a link, get the status of the source
        //      file/directory.
        //  \param maxTime The timeout value.
        //
        //  \note An empty filename is a no-op.
        fileStat
        (
            const fileName& fName,
            const bool followLink = true,
            const unsigned int maxTime = 0
        );

        //- Construct from Istream
        explicit fileStat(Istream& is);


    // Member Functions

    // Access

        //- Raw status
        const struct stat& status() const
        {
            return status_;
        }

        //- Was file-stat successful?
        bool valid() const
        {
            return valid_;
        }

        //- Size in bytes. Zero for an invalid file-stat.
        label size() const;

        //- Return the modification time in seconds.
        //  Zero for an invalid file-stat.
        time_t modTime() const;

        //- Return the modification time in seconds (nanosecond resolution)
        //  Zero for an invalid file-stat.
        double dmodTime() const;


    // Check

        //- Compare two fileStats for same device
        bool sameDevice(const fileStat& other) const;

        //- Compare two fileStats for same Inode
        bool sameINode(const fileStat& other) const;

        //- Compare state against inode
        bool sameINode(const label iNode) const;


    // IOstream Operators

        friend Istream& operator>>(Istream& is, fileStat& fs);
        friend Ostream& operator<<(Ostream& os, const fileStat& fs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
