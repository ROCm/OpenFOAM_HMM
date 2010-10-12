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

\*----------------------------------------------------------------------------*/

#include "fileMonitor.H"
#include "IOstreams.H"
#include "Pstream.H"
#include "PackedList.H"
#include "PstreamReduceOps.H"
#include "OSspecific.H"

#ifdef FOAM_USE_STAT
#   include "OSspecific.H"
#   include "regIOobject.H"     // for fileModificationSkew symbol
#else
#   include <sys/inotify.h>
#   include <sys/ioctl.h>

#   define EVENT_SIZE  ( sizeof (struct inotify_event) )
#   define EVENT_LEN   (EVENT_SIZE + 16)
#   define EVENT_BUF_LEN     ( 1024 * EVENT_LEN )
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::fileMonitor, 0);

const Foam::NamedEnum<Foam::fileMonitor::fileState, 3>
    Foam::fileMonitor::fileStateNames_;

namespace Foam
{
    template<>
    const char* Foam::NamedEnum<Foam::fileMonitor::fileState, 3>::names[] =
    {
        "unmodified",
        "modified",
        "deleted"
    };

    //- Reduction operator for PackedList of fileState
    class reduceFileStates
    {
        public:
        unsigned int operator()(const unsigned int x, const unsigned int y)
        const
        {
            // x,y are sets of 2bits representing fileState

            unsigned int mask = 3u;
            unsigned int shift = 0;
            unsigned int result = 0;

            while (mask)
            {
                // Combine state
                unsigned int xState = (x & mask) >> shift;
                unsigned int yState = (y & mask) >> shift;

                // Combine and add to result. Combine is such that UNMODIFIED
                // wins.
                unsigned int state = min(xState, yState);
                result |= (state << shift);

                shift += 2;
                mask <<= 2;
            }
            return result;
        }
    };

    //- Combine operator for PackedList of fileState
    class combineReduceFileStates
    {
        public:
        void operator()(unsigned int& x, const unsigned int y) const
        {
            x = reduceFileStates()(x, y);
        }
    };



    //! @cond internalClass
    //- Internal tracking via stat(3p) or inotify(7)
    class fileMonitorWatcher
    {
    public:

#ifdef FOAM_USE_STAT
        //- From watch descriptor to modified time
        DynamicList<time_t> lastMod_;

        //- initialize HashTable size
        inline fileMonitorWatcher(const label sz = 20)
        :
            lastMod_(sz)
        {}

        inline bool addWatch(const label watchFd, const fileName& fName)
        {
            if (watchFd < lastMod_.size() && lastMod_[watchFd] != 0)
            {
                // Reuse of watchFd : should have lastMod set to 0.
                FatalErrorIn("addWatch(const label, const fileName&)")
                    << "Problem adding watch " << watchFd
                    << " to file " << fName
                    << abort(FatalError);
            }

            lastMod_(watchFd) = lastModified(fName);
            return true;
        }

        inline bool removeWatch(const label watchFd)
        {
            lastMod_[watchFd] = 0;
            return true;
        }

#else
        //- File descriptor for the inotify instance
        int fd;

        //- Current watchIDs and corresponding directory id
        DynamicList<label> dirWatches_;
        DynamicList<fileName> dirFiles_;

        //- initialise inotify
        inline fileMonitorWatcher(const label sz = 20)
        :
            fd(inotify_init()),
            dirWatches_(sz),
            dirFiles_(sz)
        {}

        //- remove all watches
        inline ~fileMonitorWatcher()
        {
            forAll(dirWatches_, i)
            {
                if (dirWatches_[i] >= 0)
                {
                    if (inotify_rm_watch(fd, int(dirWatches_[i])))
                    {
                        WarningIn("fileMonitor::~fileMonitor()")
                            << "Failed deleting directory watch "
                            << dirWatches_[i] << endl;
                    }
                }
            }
        }

        inline bool addWatch(const label watchFd, const fileName& fName)
        {
            // Add/retrieve watch on directory containing file
            label dirWatchID = inotify_add_watch
            (
                fd,
                fName.path().c_str(),
                IN_CLOSE_WRITE
            );

            if (dirWatchID < 0)
            {
                FatalErrorIn("addWatch(const label, const fileName&)")
                    << "Failed adding watch " << watchFd
                    << " to directory " << fName
                    << exit(FatalError);
            }

            if (watchFd < dirWatches_.size() && dirWatches_[watchFd] != -1)
            {
                // Reuse of watchFd : should have dir watchID set to -1.
                FatalErrorIn("addWatch(const label, const fileName&)")
                    << "Problem adding watch " << watchFd
                    << " to file " << fName
                    << abort(FatalError);
            }

            dirWatches_(watchFd) = dirWatchID;
            dirFiles_(watchFd) = fName.name();
            return true;
        }

        inline bool removeWatch(const label watchFd)
        {
            dirWatches_[watchFd] = -1;
            return true;
        }
#endif

    };
    //! @endcond
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fileMonitor::checkFiles() const
{
#ifdef FOAM_USE_STAT
    forAll(watcher_->lastMod_, watchFd)
    {
        time_t oldTime = watcher_->lastMod_[watchFd];

        if (oldTime != 0)
        {
            const fileName& fName = watchFile_[watchFd];
            time_t newTime = lastModified(fName);

            if (newTime == 0)
            {
                state_[watchFd] = DELETED;
            }
            else
            {
                if (newTime > (oldTime + regIOobject::fileModificationSkew))
                {
                    watcher_->lastMod_[watchFd] = newTime;
                    state_[watchFd] = MODIFIED;
                }
                else
                {
                    state_[watchFd] = UNMODIFIED;
                }
            }
        }
    }
#else
    // Large buffer for lots of events
    char buffer[EVENT_BUF_LEN];

    while (true)
    {
        struct timeval zeroTimeout = {0, 0};

        //- Pre-allocated structure containing file descriptors
        fd_set fdSet;
        // Add notify descriptor to select fd_set
        FD_ZERO(&fdSet);
        FD_SET(watcher_->fd, &fdSet);

        int ready = select
        (
            watcher_->fd+1,     // num filedescriptors in fdSet
            &fdSet,             // fdSet with only inotifyFd
            NULL,               // No writefds
            NULL,               // No errorfds
            &zeroTimeout        // eNo timeout
        );

        if (ready < 0)
        {
            FatalErrorIn("fileMonitor::updateStates()")
                << "Problem in issuing select."
                << abort(FatalError);
        }
        else if (FD_ISSET(watcher_->fd, &fdSet))
        {
            // Read events
            ssize_t nBytes = read(watcher_->fd, buffer, EVENT_BUF_LEN);

            if (nBytes < 0)
            {
                FatalErrorIn("fileMonitor::updateStates(const fileName&)")
                    << "read of " << watcher_->fd
                    << " failed with " << label(nBytes)
                    << abort(FatalError);
            }

            // Go through buffer, consuming events
            int i = 0;
            while (i < nBytes)
            {
                const struct inotify_event* inotifyEvent =
                    reinterpret_cast<const struct inotify_event*>
                    (
                        &buffer[i]
                    );

                //Pout<< "watchFd:" << inotifyEvent->wd << nl
                //    << "mask:" << inotifyEvent->mask << nl
                //  << endl;
                //Pout<< "file:" << fileName(inotifyEvent->name) << endl;
                //Pout<< "len:" << inotifyEvent->len << endl;

                if ((inotifyEvent->mask & IN_CLOSE_WRITE) && inotifyEvent->len)
                {
                    // Search for file
                    forAll(watcher_->dirWatches_, i)
                    {
                        label id = watcher_->dirWatches_[i];
                        if
                        (
                            id == inotifyEvent->wd
                         && inotifyEvent->name == watcher_->dirFiles_[i]
                        )
                        {
                            // Correct directory and name
                            state_[i] = MODIFIED;
                        }
                    }
                }

                i += EVENT_SIZE + inotifyEvent->len;
            }
        }
        else
        {
            // No data
            return;
        }
    }
#endif
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::fileMonitor::fileMonitor()
:
    state_(20),
    watchFile_(20),
    freeWatchFds_(2),
    watcher_(new fileMonitorWatcher(20))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileMonitor::~fileMonitor()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::fileMonitor::addWatch(const fileName& fName)
{
    label watchFd;

    label sz = freeWatchFds_.size();
    if (sz)
    {
        watchFd = freeWatchFds_[sz-1];
        freeWatchFds_.setSize(sz-1);
    }
    else
    {
        watchFd = state_.size();
    }

    watcher_->addWatch(watchFd, fName);

    if (debug)
    {
        Pout<< "fileMonitor : added watch " << watchFd << " on file "
            << fName << endl;
    }

    if (watchFd < 0)
    {
        WarningIn("fileMonitor::addWatch(const fileName&)")
            << "could not add watch for file " << fName << endl;
    }
    else
    {
        state_(watchFd) = UNMODIFIED;
        watchFile_(watchFd) = fName;
    }
    return watchFd;
}


bool Foam::fileMonitor::removeWatch(const label watchFd)
{
    if (debug)
    {
        Pout<< "fileMonitor : removing watch " << watchFd << " on file "
            << watchFile_[watchFd] << endl;
    }

    freeWatchFds_.append(watchFd);
    return watcher_->removeWatch(watchFd);
}


const Foam::fileName& Foam::fileMonitor::getFile(const label watchFd) const
{
    return watchFile_[watchFd];
}


Foam::fileMonitor::fileState Foam::fileMonitor::getState(const label watchFd)
const
{
    return state_[watchFd];
}


void Foam::fileMonitor::updateStates(const bool syncPar) const
{
    checkFiles();

    if (syncPar)
    {
        PackedList<2> stats(state_.size());
        forAll(state_, watchFd)
        {
            stats[watchFd] = static_cast<unsigned int>(state_[watchFd]);
        }
        // Save local state for warning message below
        PackedList<2> thisProcStats(stats);

        if (stats.storage().size() == 1)
        {
            // Optimisation valid for most cases.
            reduce(stats.storage()[0], reduceFileStates());
        }
        else
        {
            Pstream::listCombineGather
            (
                stats.storage(),
                combineReduceFileStates()
            );
        }

        forAll(state_, watchFd)
        {
            if (thisProcStats[watchFd] != UNMODIFIED)
            {
                if (stats[watchFd] == UNMODIFIED)
                {
                    WarningIn("fileMonitor::updateStates(const bool) const")
                        << "Delaying reading " << watchFile_[watchFd]
                        << " due to inconsistent "
                           "file time-stamps between processors"
                        << endl;
                }
                else
                {
                    unsigned int stat = stats[watchFd];
                    state_[watchFd] = fileState(stat);
                }
            }
        }
    }
}


void Foam::fileMonitor::setUnmodified(const label watchFd)
{
#ifdef FOAM_USE_STAT
    watcher_->lastMod_[watchFd] = lastModified(watchFile_[watchFd]);
#endif
    state_[watchFd] = UNMODIFIED;
}


// ************************************************************************* //
