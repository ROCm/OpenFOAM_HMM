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

#ifdef FOAM_USE_STAT
#   include "OSspecific.H"
#   include "regIOobject.H"     // for fileModificationSkew symbol
#else
#   include <sys/inotify.h>
#   include <sys/ioctl.h>
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::fileMonitor, 0);

template<>
const char* Foam::NamedEnum<Foam::fileMonitor::fileState, 3>::names[] =
{
    "unmodified",
    "modified",
    "deleted"
};
const Foam::NamedEnum<Foam::fileMonitor::fileState, 3>
    Foam::fileMonitor::fileStateNames_;


namespace Foam
{
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
        HashTable<label, time_t> lastMod;

        //- initialize HashTable size
        inline fileMonitorWatcher(const label sz = 20)
        :
            lastMod(sz)
        {}

        inline label addWatch(const fileName& fName)
        {
            const label watchFd = lastMod.size();
            lastMod.insert(watchFd, lastModified(fName));
            return watchFd;
        }

        inline bool removeWatch(const label watchFd)
        {
            return lastMod.erase(watchFd);
        }

#else
        //- File descriptor for the inotify instance
        int fd;

        //- Pre-allocated structure containing file descriptors
        fd_set fdSet;

        //- initialize inotify
        inline fileMonitorWatcher(const label dummy = 0)
        :
            fd(inotify_init())
        {
            // Add notify descriptor to select fd_set
            FD_ZERO(&fdSet);
            FD_SET(fd, &fdSet);
        }

        //- test if file descriptor is set
        inline bool isSet() const
        {
            return FD_ISSET(fd, &fdSet);
        }

        //- reset file descriptor
        inline void reset()
        {
            FD_SET(fd, &fdSet);
        }

        inline label addWatch(const fileName& fName)
        {
            return inotify_add_watch
            (
                fd,
                fName.c_str(),
                // IN_ALL_EVENTS
                IN_CLOSE_WRITE | IN_DELETE_SELF | IN_MODIFY
            );
        }

        inline bool removeWatch(const label watchFd)
        {
            return inotify_rm_watch(fd, int(watchFd)) == 0;
        }
#endif

    };
    //! @endcond
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fileMonitor::checkFiles() const
{
#ifdef FOAM_USE_STAT
    for
    (
        HashTable<label, time_t>::iterator iter = watcher_->lastMod.begin();
        iter != watcher_->lastMod.end();
        ++iter
    )
    {
        const label watchFd = iter.key();
        const fileName& fName = watchFile_[watchFd];
        time_t newTime = lastModified(fName);

        if (newTime == 0)
        {
            state_.set(watchFd, DELETED);
        }
        else
        {
            time_t oldTime = iter();
            if (newTime > (oldTime + regIOobject::fileModificationSkew))
            {
                iter() = newTime;
                state_.set(watchFd, MODIFIED);
            }
            else
            {
                state_.set(watchFd, UNMODIFIED);
            }
        }
    }
#else
    while (true)
    {
        struct timeval zeroTimeout = {0, 0};

        int ready = select
        (
            watcher_->fd+1,        // num filedescriptors in fdSet
            &(watcher_->fdSet),    // fdSet with only inotifyFd
            NULL,                  // No writefds
            NULL,                  // No errorfds
            &zeroTimeout           // eNo timeout
        );

        if (ready < 0)
        {
            FatalErrorIn("fileMonitor::updateStates()")
                << "Problem in issuing select."
                << abort(FatalError);
        }
        else if (watcher_->isSet())
        {
            struct inotify_event inotifyEvent;

            // Read first event
            ssize_t nBytes = read
            (
                watcher_->fd,
                &inotifyEvent,
                sizeof(inotifyEvent)
            );

            if (nBytes != sizeof(inotifyEvent))
            {
                FatalErrorIn("fileMonitor::updateStates(const fileName&)")
                    << "Read " << label(nBytes) << " ; expected "
                    << label(sizeof(inotifyEvent))
                    << abort(FatalError);
            }

            // Pout<< "mask:" << inotifyEvent.mask << nl
            //     << "watchFd:" << inotifyEvent.wd << nl
            //     << "watchName:" << watchFile_[inotifyEvent.wd] << endl;

            if (inotifyEvent.mask % IN_DELETE_SELF)
            {
                Map<fileState>::iterator iter =
                    state_.find(label(inotifyEvent.wd));
                iter() = DELETED;
            }
            else if
            (
                (inotifyEvent.mask % IN_MODIFY)
             || (inotifyEvent.mask % IN_CLOSE_WRITE)
            )
            {
                Map<fileState>::iterator iter =
                    state_.find(label(inotifyEvent.wd));
                iter() = MODIFIED;
            }
        }
        else
        {
            // No data - reset
            watcher_->reset();
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
    watcher_(new fileMonitorWatcher(20))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileMonitor::~fileMonitor()
{
    // Remove watch on any remaining files
    List<label> watchFds(state_.toc());
    forAll(watchFds, i)
    {
        removeWatch(watchFds[i]);
    }

    delete watcher_;
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::fileMonitor::addWatch(const fileName& fName)
{
    const label watchFd = watcher_->addWatch(fName);

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
        state_.insert(watchFd, UNMODIFIED);
        watchFile_.insert(watchFd, fName);
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

    state_.erase(watchFd);
    watchFile_.erase(watchFd);
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
        label i = 0;
        forAllConstIter(Map<fileState>, state_, iter)
        {
            stats[i++] = static_cast<unsigned int>(iter());
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

        i = 0;
        forAllIter(Map<fileState>, state_, iter)
        {
            if (thisProcStats[i] != UNMODIFIED)
            {
                if (stats[i] == UNMODIFIED)
                {
                    WarningIn("fileMonitor::updateStates(const bool) const")
                        << "Delaying reading " << watchFile_[iter.key()]
                        << " due to inconsistent "
                           "file time-stamps between processors"
                        << endl;
                }
                else
                {
                    unsigned int stat = stats[i];
                    iter() = fileState(stat);
                }
            }
            i++;
        }
    }
}


void Foam::fileMonitor::setUnmodified(const label watchFd)
{
#ifdef FOAM_USE_STAT
    watcher_->lastMod[watchFd] = lastModified(watchFile_[watchFd]);
#endif

    Map<fileState>::iterator iter = state_.find(watchFd);

    if (iter == state_.end())
    {
        FatalErrorIn("fileMonitor::setUnmodified(const label)")
            << "Illegal watchFd " << watchFd
            << abort(FatalError);
    }

    iter() = UNMODIFIED;
}


// ************************************************************************* //
