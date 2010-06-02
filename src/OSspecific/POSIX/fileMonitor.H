/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
     \\/     M anipulation  |
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
    fileMonitor

Description
    Checking for changes to files.

!!!!!!!NOTE:
    Default is to use inotify (Linux specific, since 2.6.13)

    Compile with FOAM_USE_STAT to use the stat function call.


    - works fine except for if file gets deleted and recreated
    it stops monitoring the file!
    (does work though if the file gets moved)

SourceFiles
    fileMonitor.C

\*---------------------------------------------------------------------------*/

#ifndef fileMonitor_H
#define fileMonitor_H

#include <sys/types.h>
#include "Map.H"
#include "NamedEnum.H"
#include "labelList.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fileMonitor;

/*---------------------------------------------------------------------------*\
                           Class fileMonitor Declaration
\*---------------------------------------------------------------------------*/

class fileMonitor
{

public:

    // Public data types

        //- Enumeration defining the file state.
        enum fileState
        {
            UNMODIFIED = 0,
            DELETED = 1,
            MODIFIED = 2
        };

        static const NamedEnum<fileState, 3> fileStateNames_;

private:
    // Private data

        //- State for all watchFds
        mutable Map<fileState> state_;
        
        //- From watch descriptor to filename
        HashTable<fileName, label> watchFile_;

#ifdef FOAM_USE_STAT

        //- From watch descriptor to modified time
        mutable HashTable<label, time_t> lastModified_;
#else
        //- File descriptor for the inotify instance
        int inotifyFd_;

        //- Pre-allocated structure containing file descriptors
        mutable fd_set watchSet_;

#endif
        //- Update state_ from any events.
        void checkFiles() const;

        //- Disallow default bitwise copy construct
        fileMonitor(const fileMonitor&);

        //- Disallow default bitwise assignment
        void operator=(const fileMonitor&);


public:

    // Declare name of the class and its debug switch
    ClassName("fileMonitor");

    // Constructors

        //- Construct null
        fileMonitor();


    // Destructor

        ~fileMonitor();


    // Member Functions

        //- Add file to watch. Returns watch descriptor
        label addWatch(const fileName&);

        //- Remove file to watch. Return true if successful
        bool removeWatch(const label watchFd);

        //- Get name of file being watched
        const fileName& getFile(const label watchFd) const;

        //- Check state using handle
        fileState getState(const label watchFd) const;

        //- Check state of all files. Updates state_.
        void updateStates(const bool syncPar) const;

        //- Reset state (e.g. after having read it) using handle
        void setUnmodified(const label watchFd);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
